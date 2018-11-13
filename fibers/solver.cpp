#include <set>
#include <valarray>

#include "solver.h"
#include "constants.h"
#include "option.h"

void getVert2Elements(Mesh *mesh, vector<vector<int> >& vert2Elements) {

    int NumOfVertices = mesh->GetNV();
    for (int i = 0; i < NumOfVertices; i++) {
        vector<int> elements;
        vert2Elements.push_back(elements);
    }

    int NumOfElements = mesh->GetNE();
    for (int i = 0; i < NumOfElements; i++) {
        const Element *ele = mesh->GetElement(i);
        const int *v = ele->GetVertices();
        const int nv = ele->GetNVertices();
        for (int j = 0; j < nv; j++) {
            int vert = v[j];
            vert2Elements[vert].push_back(i);
        }
    }
    
    for (int i = 0; i < NumOfVertices; i++) {
        stringstream msg;
        msg << "getVert2Elements : vertex[" << i << "] size is zero"; 
        MFEM_ASSERT(vert2Elements[i].size()!=0, msg.str());
    }
    

}


GridFunction laplace(Mesh *mesh, Array<int> &all_ess_bdr, Array<int> &nonzero_ess_bdr, Array<int> &zero_ess_bdr, Option& options, int myid){
    
    int dim = mesh->Dimension();
    if(options.verbose){
        cout << "\tDimension =" << dim << endl;    
    }
    // 4. Define a finite element space on the mesh. Here we use continuous
    //    Lagrange finite elements of the specified order. If order < 1, we
    //    instead use an isoparametric/isogeometric space.
    FiniteElementCollection *fec;
    if (options.order > 0) {
        fec = new H1_FECollection(options.order, dim);
    } else if (mesh->GetNodes()) {
        fec = mesh->GetNodes()->OwnFEC();
        cout << "\tUsing isoparametric FEs: " << fec->Name() << endl;
    } else {
        fec = new H1_FECollection(options.order = 1, dim);
    }
    FiniteElementSpace *fespace = new FiniteElementSpace(mesh, fec);
    if (myid == 0 && options.verbose) {
        cout << "\tNumber of finite element unknowns: "
                << fespace->GetTrueVSize() << endl;
    }

    // 5. Determine the list of true (i.e. conforming) essential boundary dofs.
    //    In this example, the boundary conditions are defined by marking all
    //    the boundary attributes from the mesh as essential (Dirichlet) and
    //    converting them to a list of true dofs.
    Array<int> ess_tdof_list;
    MFEM_ASSERT(mesh->bdr_attributes.Size()!=0, "Boundary size cannot be zero."); 
    
    fespace->GetEssentialTrueDofs(all_ess_bdr, ess_tdof_list);
  

    // 6. Set up the linear form b(.) which corresponds to the right-hand side of
    //    the FEM linear system, which in this case is (1,phi_i) where phi_i are
    //    the basis functions in the finite element fespace.
    LinearForm *b = new LinearForm(fespace);
    ConstantCoefficient zero(0.0);
    b->AddDomainIntegrator(new DomainLFIntegrator(zero));
    b->Assemble();

    // 7. Define the solution vector x as a finite element grid function
    //    corresponding to fespace. Initialize x with initial guess of zero,
    //    which satisfies the boundary conditions.
    GridFunction x(fespace);
    x = 0.0;

    if (myid == 0 && options.verbose) {
        cout << "\tall_ess_bdr size=" << all_ess_bdr.Size() << endl;
        cout << "\tx size " << x.Size() << endl;
    }

    // 8. Set up the bilinear form a(.,.) on the finite element space
    //    corresponding to the Laplacian operator -Delta, by adding the Diffusion
    //    domain integrator.
    BilinearForm *a = new BilinearForm(fespace);

    // The diffusion integrator should have a coefficient of one, not zero
    ConstantCoefficient one(1.0);
    a->AddDomainIntegrator(new DiffusionIntegrator(one));

    // 9. Assemble the bilinear form and the corresponding linear system,
    //    applying any necessary transformations such as: eliminating boundary
    //    conditions, applying conforming constraints for non-conforming AMR,
    //    static condensation, etc.
    if (options.static_cond) {
        a->EnableStaticCondensation();
    }
    a->Assemble();

    // Project the constant 14 value to all boundary attributes except 1
    ConstantCoefficient nonzero_bdr(1.0);
    x.ProjectBdrCoefficient(nonzero_bdr, nonzero_ess_bdr);

    // Project the constant 0 value to boundary attribute 1
    ConstantCoefficient zero_bdr(0.0);
    x.ProjectBdrCoefficient(zero_bdr, zero_ess_bdr);
        
    //return 0;
    SparseMatrix A;
    Vector B, X;
    // Form the linear system using ALL of the essential boundary dofs (from all_ess_bdr)
    a->FormLinearSystem(ess_tdof_list, x, *b, A, X, B);

    if (myid == 0 && options.verbose) {
        cout << "\tSize of linear system: " << A.Height() << endl;
    }
#ifndef MFEM_USE_SUITESPARSE
    // 10. Define a simple symmetric Gauss-Seidel preconditioner and use it to //ILU
    //     solve the system A X = B with PCG.
    GSSmoother M(A);
    if (myid == 0) {
        int printLevel=-1;
        if(options.verbose){
              printLevel=1;  
        }            
        PCG(A, M, B, X, printLevel, 1000, 1e-12, 0.0);
    }else{
        // always turn of print out for rank>0
        PCG(A, M, B, X, -1, 1000, 1e-12, 0.0);
    }
#else
    // 10. If MFEM was compiled with SuiteSparse, use UMFPACK to solve the system.
    UMFPackSolver umf_solver;
    umf_solver.Control[UMFPACK_ORDERING] = UMFPACK_ORDERING_METIS;
    umf_solver.SetOperator(A);
    umf_solver.Mult(B, X);
#endif

    // 11. Recover the solution as a finite element grid function.
    a->RecoverFEMSolution(X, *b, x);
      
    return x;
}


void getVetecesGradients(Mesh *mesh, GridFunction& x, vector<vector<int> >& vert2Elements, vector<double> &pot, vector<Vector> &gradients, string output, int myid){
    //double *x_data=x.GetData();
    for(int i=0; i<x.Size(); i++){         
        //pot.push_back(x_data[i]);
        pot.push_back(x(i));
        //if(i<5) cout << "pot " << x(i) << endl;
    }
    
    const FiniteElementSpace *fes=x.FESpace();
    unsigned nv = fes->GetNV();
    for (unsigned i = 0; i < nv; i++) {
        vector<int> elements = vert2Elements[i];
        Vector grad(3);
        Vector grad_ele(3);
        grad=0.0;
        stringstream msg;
        msg << "laplace : vertex[" << i << "] size is zero";
        MFEM_ASSERT(elements.size()!=0, msg.str());
        for (unsigned j = 0; j < elements.size(); j++) {
            ElementTransformation * tr = fes->GetElementTransformation(elements[j]);
            const IntegrationRule &ir = fes->GetFE(elements[j])->GetNodes();  // Get the parametric integration rule
            grad_ele= 0.0;
            for (int k=0; k < ir.GetNPoints(); k++) {
               Vector grad_point(3);
               grad_point = 0.0;
               const IntegrationPoint &ip = ir.IntPoint(k); // Get the current integration point
               tr->SetIntPoint(&ip); // Set the integration point for the transformation
               x.GetGradient((*tr), grad_point);
               grad_ele += grad_point;
            }
            grad_ele /= ir.GetNPoints();
            grad+=grad_ele;
        }
        grad/=elements.size();
        gradients.push_back(grad);

//        if (i < 5) {
//            cout << "grad ";
//            for (int j = 0; j < grad.Size(); j++) {
//                cout << grad[j] << " ";
//            }
//            cout << endl;
//        }

    }

    // 12. Save the refined mesh and the solution. This output can be viewed later
    //     using GLVis: "glvis -m refined.mesh -g sol.gf".

    if (myid == 0) {
        string fileName = output + ".mesh";
        ofstream mesh_ofs(fileName.c_str());
        mesh_ofs.precision(8);
        mesh->Print(mesh_ofs);

        fileName = output + ".gf";
        ofstream sol_ofs(fileName.c_str());
        sol_ofs.precision(8);
        x.Save(sol_ofs);

        vector<Vector> xvectors;
        for (int i = 0; i < x.Size(); i++) {
            Vector vectmp(3);
            vectmp(0) = x(i);
            vectmp(1) = 0.0;
            vectmp(2) = 0.0;
            xvectors.push_back(vectmp);
        }

        fileName = output + "-x.vtk";
        ofstream x_ofs(fileName.c_str());
        printFiberVTK(mesh, xvectors, x_ofs);
    }

    // 14. Free the used memory.
//    delete a;
//    delete b;
//    delete fespace;
//    if (order > 0) {
//        delete fec;
//    }
    
}

bool isPlanar(double *coor0, double *coor1, double *coor2, double cosTheta, int axis){
    double u[3];
    double v[3];
    double w[3];
    u[0]=coor1[0]-coor0[0];
    u[1]=coor1[1]-coor0[1];
    u[2]=coor1[2]-coor0[2];
    v[0]=coor2[0]-coor0[0];
    v[1]=coor2[1]-coor0[1];
    v[2]=coor2[2]-coor0[2];
    
    w[0]=u[1]*v[2]-u[2]*v[1];
    w[1]=u[2]*v[0]-u[0]*v[2];
    w[2]=u[0]*v[1]-u[1]*v[0];
    
    double r=sqrt(w[0]*w[0]+w[1]*w[1]+w[2]*w[2]);
    double cosT=abs(w[axis]/r);
    
    if(cosT>cosTheta){
        return true;
    }
    
    return false;
    
}

bool isTriInTet(vector<int>& tri, vector<int>& tet){
    MFEM_ASSERT(tri.size() == 3, "Wrong boundary size");
    MFEM_ASSERT(tet.size() == 4, "Wrong tetrahedral size");

    if(tri[0]==tet[0]){
        if(tri[1]==tet[1] && tri[2]==tet[2]){
            return true;
        }        
    }else if(tri[0]==tet[1]){
        if(tri[1]==tet[2] && tri[2]==tet[3]){
            return true;
        }          
    }    
    return false;
}

void findNeighborRecursive(Element* ele, vector<Element*>& elements, int attr){
    const int *v = ele->GetVertices();
    const int nv = ele->GetNVertices(); 
    for(unsigned i=0; i<elements.size(); i++){
        Element* queryEle=elements[i];
        // Only search for elements with unassigned attributes.
        if(queryEle->GetAttribute()==0){ 
            const int *qv = queryEle->GetVertices();
            const int nqv = queryEle->GetNVertices(); 
            bool isNeighbor=false;
            for (int j = 0; j < nv; j++) {
                for (int k = 0; k < nqv; k++) {
                    // If two elements share the same vertex they are neighbor. 
                    if(v[j]==qv[k]){                   
                        isNeighbor=true;
                        // Should break two loops can use lambda or function return.
                        break;  
                    }
                }
            }
            if(isNeighbor){
                queryEle->SetAttribute(attr);
                // recursively search for neighboring elements.
                findNeighborRecursive(queryEle, elements, attr);
            }
        }
    }           
}

set<Element*> findNeighbor(Element* ele, set<Element*>& elements, int attr){
    
    set<Element*> explored;
    set<Element*> front;
    
    front.insert(ele);
    
    while(!front.empty()){
        // grab item from front
        Element* item=*(front.begin());        
        explored.insert(item);
        front.erase(front.begin());
        
        set<Element*> neighbors;

        const int *v = item->GetVertices();
        const int nv = item->GetNVertices(); 
        for(set<Element*>::iterator it=elements.begin(); it!=elements.end(); ++it){
            Element* queryEle=(*it);
            const int *qv = queryEle->GetVertices();
            const int nqv = queryEle->GetNVertices(); 
            bool found=false;
            for (int j = 0; j < nv; j++) {
                for (int k = 0; k < nqv; k++) {
                    // If two elements share the same vertex they are neighbor. 
                    if(v[j]==qv[k]){                   
                        neighbors.insert(queryEle);
                        found=true;
                        //cout << "DEBUG size of neighbors=" << neighbors.size()<< endl;
                        break;
                    }
                }
                if(found) break;
            }
            
        }
                        
        for(set<Element*>::iterator it=neighbors.begin(); it!=neighbors.end(); ++it){
            elements.erase(*(it)); // remove the found neighbors from elements
            
            if(explored.find(*(it)) == explored.end()){ // if neighobr is not in explored
                front.insert(*(it));
            }
        }
               
    }
    
    for(set<Element*>::iterator it=explored.begin(); it!=explored.end(); ++it){
        (*it)->SetAttribute(attr); //set the attribute values
    }
    
    //cout << "DEBUG number of explored elements=" << explored.size() << endl;
    
    return explored;
}

int detectAxis(Mesh *mesh, vector<Vector>& boundingbox){
    // 0 - x; 1 - y; 2 - z
    assert(boundingbox.size()==2);
    Vector coord_min=boundingbox[0];
    Vector coord_max=boundingbox[1];
    
    vector<double> top5(3);
    
    for(int i=0; i<3; i++){
        top5[i]=coord_max[i]-(coord_max[i]-coord_min[i])*0.05;
    }
    
    vector<int> count(3, 0);
    
    int nbe=mesh->GetNBE();
    for(int i=0; i<nbe; i++){
        Element *ele = mesh->GetBdrElement(i);        
        const int *v = ele->GetVertices();
        const int nv = ele->GetNVertices();
        MFEM_ASSERT(nv == 3, "Wrong boundary size");
        
        double *coord0 = mesh->GetVertex(v[0]);
        
        for(int i=0; i<3; i++){
            if(coord0[i]>top5[i]){
                count[i]++;
            }            
        }

    }

    int axis=0;
    int maxVal=-1;
    for(int i=0; i<3; i++){
        if(count[i]>maxVal){
            maxVal=count[i];
            axis=i;
        }
    }
    
    return axis; 
    
}

void setSurfaces(Mesh *mesh, vector<Vector>& boundingbox, double angle, int myid){
    // Attributes for different surface
    const int apexAttr=1;
    const int baseAttr=2;
    const int epiAttr=3;
    const int lvAttr=5; 
    const int rvAttr=4;
       
    // Determine the max and min dimension of mesh and apex.
    double *coord;
    Vector coord_min(3);
    Vector coord_max(3);
    bool firstEle=true;
    
    int nbe=mesh->GetNBE();
    for(int i=0; i<nbe; i++){
        Element *ele = mesh->GetBdrElement(i);        
        const int *v = ele->GetVertices();
        const int nv = ele->GetNVertices();
        // The first loop has to initialize the min and max.
        if(firstEle){
            firstEle=false;
            coord=mesh->GetVertex(v[0]);
            for (int j = 0; j < 3; j++) {
                coord_min(j)=coord[j];
                coord_max(j)=coord[j];
            }            
        }
        
        for(int j=0; j<nv; j++){
            coord=mesh->GetVertex(v[j]);
            
            for (int k = 0; k < 3; k++) {
                if(coord[k]<coord_min[k]){
                    coord_min(k)=coord[k];
                }
                if(coord[k]>coord_max[k]){
                    coord_max(k)=coord[k];
                }            
            }                                    
        }
        
    }
    
    boundingbox.push_back(coord_min);
    boundingbox.push_back(coord_max);
    
    // axis: 0 - x; 1 - y; 2 - z
    int axis=detectAxis(mesh, boundingbox);

    // Keep track vertex and element indeces for min in axis
    int apexVet=0;
    int apexEleIndex=0;    
    double apexCoord_min=0;
    firstEle=true;
    
    for(int i=0; i<nbe; i++){
        Element *ele = mesh->GetBdrElement(i);        
        const int *v = ele->GetVertices();
        const int nv = ele->GetNVertices();
        
        // The first loop has to initialize the min and max.
        if(firstEle){
            firstEle=false;
            coord=mesh->GetVertex(v[0]);
            apexCoord_min=coord[axis];
            apexVet = v[0];
            apexEleIndex = i;
        }
        
        for(int j=0; j<nv; j++){
            coord=mesh->GetVertex(v[j]);
                       
            if(coord[axis]<apexCoord_min){
                apexCoord_min=coord[axis];  
                apexVet=v[j];
                apexEleIndex=i;
                
            }                                   
        }
        
    }    

    coord = mesh->GetVertex(apexVet);
    
    // Top 5% of the  axis.
    double apexTop5=coord_max[axis]-(coord_max[axis]-coord_min[axis])*0.05;

    if (myid == 0) {
        cout << "\tHeart aligns along " << axis << " (0=x, 1=y, 2=z)" << endl;
        cout << "\tMin: " << coord_min(0) << " " << coord_min(1) << " " << coord_min(2) << endl;
        cout << "\tMax: " << coord_max(0) << " " << coord_max(1) << " " << coord_max(2) << endl;
        cout << "\tApex: " << coord[0] << " " << coord[1] << " " << coord[2] << endl;
        cout << "\tTop 5% axis {"<< axis <<") coordinate: " << apexTop5 << endl;
    }
    // Initialization the attributes to 0 and set attribute of apex
    for(int i=0; i<nbe; i++){
        Element *ele = mesh->GetBdrElement(i);        
        const int *v = ele->GetVertices();
        const int nv = ele->GetNVertices();
        // initialize the attribute for boundary.  
        ele->SetAttribute(0);
        
        //Found apex elements and set attribute.
        for (int j = 0; j < nv; j++) {
            if (v[j] ==apexVet){
                ele->SetAttribute(apexAttr);
                //cout << "Element index = " << i << endl;
            }
        }        
    }
    
    // Base    
    // The base must be planar. Its norm must be within 20 degrees of z axis.
    double cosTheta = cos(angle*PI/180); 
    for(int i=0; i<nbe; i++){
        Element *ele = mesh->GetBdrElement(i);        
        const int *v = ele->GetVertices();
        const int nv = ele->GetNVertices();
        MFEM_ASSERT(nv == 3, "Wrong boundary size");
        
        double *coord0 = mesh->GetVertex(v[0]);
        if(coord0[axis]>apexTop5){
            double *coord1 = mesh->GetVertex(v[1]);
            double *coord2 = mesh->GetVertex(v[2]);
            if(isPlanar(coord0, coord1, coord2, cosTheta, axis)){
                ele->SetAttribute(baseAttr);
            }
        }
    }
    
    //EPI
    set<Element*> elements;
    for(int i=0; i<nbe; i++){
        Element *ele = mesh->GetBdrElement(i); 
        if(ele->GetAttribute()==0){
            elements.insert(ele);
        }
    }
    
    Element *apexEle=mesh->GetBdrElement(apexEleIndex);
    set<Element*> epiElements=findNeighbor(apexEle, elements, epiAttr);
        
    // set the apexEle back to apexAttr, it got changed in the findNeighbor
    apexEle->SetAttribute(apexAttr);
    
    // LV & RV
    // pick one element in the container and assume it is right ventricles.
    Element* lastEle=*(elements.begin()); 
    lastEle->SetAttribute(rvAttr);
    // get rid of last element in the container
    elements.erase(elements.begin());
    set<Element*> vElements=findNeighbor(lastEle, elements, rvAttr);
        
    //get the numbers of points in two ventricles
    unsigned v1count=vElements.size(); // assume right ventricle
    unsigned v2count=elements.size();  // assume left ventricle
 
    if (myid == 0) {
        cout << "\nEPI number of element = " << epiElements.size() << endl;
        cout << "V1  number of element = " << v1count << endl;
        cout << "V2  number of element = " << v2count << endl << endl;
    }

    //The right ventricle has more points/cells than the left 
    if(v1count>v2count){ //the assumption validated
        for(set<Element*>::iterator it=elements.begin(); it!=elements.end(); ++it){
            (*it)->SetAttribute(lvAttr); //set the attribute values
        }              
    }else{
       lastEle->SetAttribute(lvAttr);
        for(set<Element*>::iterator it=vElements.begin(); it!=vElements.end(); ++it){
            (*it)->SetAttribute(lvAttr); //set the attribute values
        }
        for(set<Element*>::iterator it=elements.begin(); it!=elements.end(); ++it){
            (*it)->SetAttribute(rvAttr); //set the attribute values
        }              
    }

    // Check if there are unassigned elements left.
    for(int i=0; i<nbe; i++){
        Element *ele = mesh->GetBdrElement(i); 
        MFEM_ASSERT(ele->GetAttribute()!=0, "Unassigned element.");
    }  
    
    mesh->SetAttributes();
            
}

void setSurf4Surf(Mesh *surface, double angle){
    // Attributes for different surface
    const int apexAttr=1;
    const int baseAttr=2;
    const int epiAttr=3;
    const int lvAttr=5; 
    const int rvAttr=4;
       
    // Determine the max and min dimension of mesh and apex.
    double *coord;
    Vector coord_min(3);
    Vector coord_max(3);
    bool firstEle=true;
    
    int ne=surface->GetNE();
    for(int i=0; i<ne; i++){
        Element *ele = surface->GetElement(i);        
        const int *v = ele->GetVertices();
        const int nv = ele->GetNVertices();
        // The first loop has to initialize the min and max.
        if(firstEle){
            firstEle=false;
            coord=surface->GetVertex(v[0]);
            for (int j = 0; j < 3; j++) {
                coord_min(j)=coord[j];
                coord_max(j)=coord[j];
            }            
        }
        
        for(int j=0; j<nv; j++){
            coord=surface->GetVertex(v[j]);
            
            for (int k = 0; k < 3; k++) {
                if(coord[k]<coord_min[k]){
                    coord_min(k)=coord[k];
                }
                if(coord[k]>coord_max[k]){
                    coord_max(k)=coord[k];
                }            
            }                                    
        }
        
    }
    
    vector<Vector> boundingbox;
    boundingbox.push_back(coord_min);
    boundingbox.push_back(coord_max);
    
    // axis: 0 - x; 1 - y; 2 - z
    int axis=detectAxis(surface, boundingbox);

    // Keep track vertex and element indeces for min in axis
    int apexVet=0;
    int apexEleIndex=0;    
    double apexCoord_min=0;
    firstEle=true;
    
    for(int i=0; i<ne; i++){
        Element *ele = surface->GetBdrElement(i);        
        const int *v = ele->GetVertices();
        const int nv = ele->GetNVertices();
        
        // The first loop has to initialize the min and max.
        if(firstEle){
            firstEle=false;
            coord=surface->GetVertex(v[0]);
            apexCoord_min=coord[axis];
            apexVet = v[0];
            apexEleIndex = i;
        }
        
        for(int j=0; j<nv; j++){
            coord=surface->GetVertex(v[j]);
                       
            if(coord[axis]<apexCoord_min){
                apexCoord_min=coord[axis];  
                apexVet=v[j];
                apexEleIndex=i;
                
            }                                   
        }
        
    }   
    
    coord = surface->GetVertex(apexVet);
    

    // Top 5% of the  axis.
    double apexTop5=coord_max[axis]-(coord_max[axis]-coord_min[axis])*0.05;


    cout << "\tHeart aligns along " << axis << " (0=x, 1=y, 2=z)" << endl;
    cout << "\tMin: " << coord_min(0) << " " << coord_min(1) << " " << coord_min(2) << endl;
    cout << "\tMax: " << coord_max(0) << " " << coord_max(1) << " " << coord_max(2) << endl;
    cout << "\tApex: " << coord[0] << " " << coord[1] << " " << coord[2] << endl;
    cout << "\tTop 5% axis {"<< axis <<") coordinate: " << apexTop5 << endl;
       
    
    // Initialization the attributes to 0 and set attribute of apex
    for(int i=0; i<ne; i++){
        Element *ele = surface->GetElement(i);        
        const int *v = ele->GetVertices();
        const int nv = ele->GetNVertices();
        // initialize the attribute for boundary.  
        ele->SetAttribute(0);
        
        //Found apex elements and set attribute.
        for (int j = 0; j < nv; j++) {
            if (v[j] ==apexVet){
                ele->SetAttribute(apexAttr);
                //cout << "Element index = " << i << endl;
            }
        }        
    }
    
    // Base    
    // The base must be planar. Its norm must be within 20 degrees of z axis.
    double cosTheta = cos(angle*PI/180); 
    for(int i=0; i<ne; i++){
        Element *ele = surface->GetElement(i);        
        const int *v = ele->GetVertices();
        const int nv = ele->GetNVertices();
        MFEM_ASSERT(nv == 3, "Wrong boundary size");
        
        double *coord0 = surface->GetVertex(v[0]);
        if(coord0[axis]>apexTop5){
            double *coord1 = surface->GetVertex(v[1]);
            double *coord2 = surface->GetVertex(v[2]);
            if(isPlanar(coord0, coord1, coord2, cosTheta, axis)){
                ele->SetAttribute(baseAttr);
            }
        }
    }
    
    //EPI
    set<Element*> elements;
    for(int i=0; i<ne; i++){
        Element *ele = surface->GetElement(i); 
        if(ele->GetAttribute()==0){
            elements.insert(ele);
        }
    }
    
    Element *apexEle=surface->GetElement(apexEleIndex);
    set<Element*> epiElements=findNeighbor(apexEle, elements, epiAttr);
    
    // set the apexEle back to apexAttr, it got changed in the findNeighbor
    apexEle->SetAttribute(apexAttr);    
    
    // LV & RV
    // pick one element in the container and assume it is right ventricles.
    Element* lastEle=*(elements.begin()); 
    lastEle->SetAttribute(rvAttr);
    // get rid of last element in the container
    elements.erase(elements.begin());
    set<Element*> vElements=findNeighbor(lastEle, elements, rvAttr);
    
    //get the numbers of points in two ventricles
    unsigned v1count=vElements.size(); // assume right ventricle
    unsigned v2count=elements.size();  // assume left ventricle
    
    //The right ventricle has more points/cells than the left 
    if(v1count>v2count){ //the assumption validated
        for(set<Element*>::iterator it=elements.begin(); it!=elements.end(); ++it){
            (*it)->SetAttribute(lvAttr); //set the attribute values
        }              
    }else{
       lastEle->SetAttribute(lvAttr);
        for(set<Element*>::iterator it=vElements.begin(); it!=vElements.end(); ++it){
            (*it)->SetAttribute(lvAttr); //set the attribute values
        }
        for(set<Element*>::iterator it=elements.begin(); it!=elements.end(); ++it){
            (*it)->SetAttribute(rvAttr); //set the attribute values
        }              
    }

    // Check if there are unassigned elements left.
    for(int i=0; i<ne; i++){
        Element *ele = surface->GetElement(i); 
        MFEM_ASSERT(ele->GetAttribute()!=0, "Unassigned element.");
    }  
    
    //surface->SetAttributes();
            
}

/*
void setBaseOLD(Mesh *mesh, int attr){
    int nv = mesh->GetNV();
    double *coord;
    double coord_min[3];
    double coord_max[3];
    int apexVet=0;
    coord = mesh->GetVertex(0);
    if(coord!=NULL){
        for (int j = 0; j < 3; j++) {
            coord_min[j]=coord[j];
            coord_max[j]=coord[j];
        }
    }
    for (int i = 0; i < nv; i++) {
        coord = mesh->GetVertex(i);
        for (int j = 0; j < 3; j++) {
            if(coord[j]<coord_min[j]){
                coord_min[j]=coord[j];
                if(j==2){
                    apexVet=i;
                }
            }
            if(coord[j]>coord_max[j]){
                coord_max[j]=coord[j];
            }            
        }

    }

    cout << "Min: " << coord_min[0] << " " << coord_min[1] << " " << coord_min[2] << endl;
    cout << "Max: " << coord_max[0] << " " << coord_max[1] << " " << coord_max[2] << endl;
    coord = mesh->GetVertex(apexVet);
    cout << "Apex: " << coord[0] << " " << coord[1] << " " << coord[2] << endl;
    
    // Top 5% of the z axis.
    double zTop5=coord_max[2]-(coord_max[2]-coord_min[2])*0.05;
    cout << "Top 5% z coordinate: " << zTop5 << endl;
    
    
    const int apexAttr=0;
    int ne = mesh->GetNE();    
    for (int i = 0; i < ne; i++) {
        Element *ele = mesh->GetElement(i);
        ele->SetAttribute(5);
        const int *v = ele->GetVertices();
        const int nv = ele->GetNVertices();
        for (int j = 0; j < nv; j++) {
            if (v[j] ==apexVet){
                ele->SetAttribute(apexAttr);
                cout << "Element index = " << i << endl;
            }
        }
    }    
    
    // for debug
    ofstream befh;
    befh.open("boundary3.txt");
    //
    
    double cosTheta = cos(20*PI/180); // within 20 degrees of z axis.
    vector<vector<int> > baseBoundary;
    int nbe=mesh->GetNBE();
    for(int i=0; i<nbe; i++){
        Element *ele = mesh->GetBdrElement(i);        
        const int *v = ele->GetVertices();
        const int nv = ele->GetNVertices();
        if(nv!=3){
            cout << "Boundary element should be 3 but it is " << nv <<endl;
            return;
        }
        double *coord0 = mesh->GetVertex(v[0]);
        if(coord0[2]>zTop5){
            double *coord1 = mesh->GetVertex(v[1]);
            double *coord2 = mesh->GetVertex(v[2]);
            if(isPlanar(coord0, coord1, coord2, cosTheta)){
                //baseBoundary.Append(ele);
                vector<int> vertecies;                
                for (int j = 0; j < nv; j++) {
                    vertecies.push_back(v[j]);
                }
                sort(vertecies.begin(), vertecies.end());
                baseBoundary.push_back(vertecies);
                // for debug
                befh << nv;
                for (int j = 0; j < nv; j++) {                    
                    befh << ' ' << v[j];                    
                }
                befh << '\n'; 
                
                befh << nv;
                for (int j = 0; j < nv; j++) {                    
                    befh << ' ' << vertecies[j];                    
                }
                befh << '\n';                 
                //
 
            }
        }
        
    }
      
    for (int i = 0; i < ne; i++) {
        Element *ele = mesh->GetElement(i);        
        const int *v = ele->GetVertices();
        const int nv = ele->GetNVertices();
        double *coord0 = mesh->GetVertex(v[0]);        
        if(coord0[2]>zTop5){
            vector<int> tet;
            for (int j = 0; j < nv; j++) {
                tet.push_back(v[j]);
            }
            sort(tet.begin(), tet.end());

            for(unsigned j=0; j < baseBoundary.size(); j++){
                vector<int> tri=baseBoundary[j];
                if(isTriInTet(tri, tet)){
                    ele->SetAttribute(attr);
                    break;
                }
            }
        }
        
    }    
    
    
    
}
*/
