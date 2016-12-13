function [thisTolerance, hornNumerator, hornDenominator] = rationalApprox(fff, outputName, inputName, lb, ub, numPoints, tolerance, errorType)
  useHornApprox = 1;
  zzz = linspace(-1,1,numPoints)';
  xxx = zzz*(ub-lb)/2 + (ub+lb)/2;
  yyy = fff(xxx);

  for cost=3:64
    [numerator, denominator, thisTolerance] = findBestZAtCost(yyy, zzz, cost, errorType);
    numTermCount = rows(numerator);
    denomTermCount = rows(denominator);
    thisTolerance

    if (thisTolerance < tolerance)
      %print it out
      if (useHornApprox)
        [numerator, denominator] = findHornFromZ(numerator,denominator,lb,ub);
      endif
      
      numeratorString = cell(size(numerator));
      denominatorString = cell(size(denominator));
      for ii=1:length(numerator)
        numeratorString{ii} = sprintf("%.17e", numerator(ii));
      endfor
      for ii=1:length(denominator)
        denominatorString{ii} = sprintf("%.17e", denominator(ii));
      endfor
      printf("double %s;\n", outputName);
      printf("{\n");
      printf("   //rationalApprox(%s, \"%s\", \"%s\", %f, %f, %d, %e, \"%s\")\n", func2str(fff), outputName, inputName, lb, ub, numPoints, tolerance, errorType);
      if (useHornApprox)
        printf("   double __x = %s;\n", inputName);
      else
        printf("   double __x = (2/(%f - %f))*(%s) - (%f + %f)/(%f - %f);\n",ub,lb,inputName,ub,lb,ub,lb);
      endif
      printf("   if (useRatpolyApprox && %f <= __x && __x <= %f)\n", lb, ub);
      printf("   {\n");
      printf("      const int pSize = %d;\n", numTermCount);
      printf("      const double pCoeff[] = { %s }\;\n", strjoin(numeratorString, ", "));
      printf("      double numerator=pCoeff[pSize-1];\n");
      printf("      for (int ip=pSize-2; ip>=0; --ip)\n");
      printf("      {\n");
      printf("         numerator = pCoeff[ip] + __x*numerator;\n")
      printf("      }\n");
      if (denomTermCount == 1)
        printf("      %s = numerator;\n", outputName);
      else
        printf("      const int qSize = %d;\n", denomTermCount);
        printf("      double qCoeff[] = { %s };\n", strjoin(denominatorString, ", "));
        printf("      double denominator = qCoeff[qSize-1];\n");
        printf("      for (int iq=qSize-2; iq>=0; --iq)\n");
        printf("      {\n");
        printf("         denominator = qCoeff[iq] + __x*denominator;\n");
        printf("      }\n");
        printf("      %s = numerator/denominator;\n", outputName);
      endif
      printf("   }\n");
      printf("   else\n");
      printf("   {\n\n");
      printf("      %s = ;\n", outputName);
      printf("   }\n");
      printf("}\n");
      break;
    endif
  endfor
endfunction

function [zNumeratorCoeffs, zDenominatorCoeffs, bestTolerance] = findBestZAtCost(yyy, zzz, maxTermCount, errorType)
  % build the polynomial
  assert(maxTermCount >= 2)
  chebyPoly = zeros(rows(zzz), maxTermCount);
  chebyPoly(:,1) = 1;
  chebyPoly(:,2) = zzz;
  for ii=3:maxTermCount
    chebyPoly(:,ii) = 2.*zzz.*chebyPoly(:,ii-1) - chebyPoly(:,ii-2);
  endfor
  figure(1);
  plot(zzz, chebyPoly);

  funcMin = min(yyy);
  funcMax = max(yyy);
  
  bestChebyCoeffs = zeros(maxTermCount-1);
  bestTolerance = 1e30;
  bestNumTermCount = 0;
  bestDenomTermCount = 0;
  for denomTermCount=1:maxTermCount-1
    numTermCount = maxTermCount - denomTermCount;
    
    [approximation, chebyCoeffs] = buildApproximate(yyy, chebyPoly, numTermCount, denomTermCount);
    rmsErrorScaled = norm(approximation - yyy)/sqrt(rows(yyy))/(funcMax-funcMin);
    maxErrorScaled = max(abs(approximation - yyy))/(funcMax-funcMin);
    if (errorType == "rms")
      errorScaled = rmsErrorScaled;
    else
      errorScaled = maxErrorScaled;
    endif
    if (errorScaled < bestTolerance)
      bestTolerance = errorScaled;
      bestNumTermCount = numTermCount;
      bestDenomTermCount = denomTermCount;
      bestChebyCoeffs = chebyCoeffs;
      figure(3)
      plot(zzz, [yyy approximation]);
      figure(4)
      plot(zzz, approximation-yyy);
    endif

    figure(2);
    plot(zzz, [yyy approximation]);
  endfor

  % Comments are in order here.
  % For chebychev, we have
  % c0 =  z0            = 1
  % c1 =      z1        = z
  % c2 = -z0     +2*z2  = 2z^2 - 1
  %
  % For horn, we have
  % x0 =   z0                     = 1
  % x1 = b*z0   +   m*z1          = b+mz
  % x2 = b^2*z0 + 2mb*z1 + m^2*z2 = b^2+2mbz+m^2*z^2
  %
  % we just found coeffs q such that q'*C is a good approx.
  % From the above, we know
  % C = DZ
  % X = BZ
  % Therefore
  % q'*C = q'*D*Z
  %      = q'*D*inv(B)*X
  %      = (inv(B')*D'*q)'*X
  
  
  chebyPolyFromZ = zeros(maxTermCount,maxTermCount);
  chebyPolyFromZ(1:2,1:2) = eye(2);
  for ii=3:maxTermCount
    chebyPolyFromZ(ii,:) = [0 2*chebyPolyFromZ(ii-1, 1:end-1)] - chebyPolyFromZ(ii-2,:);
  endfor

  zCoeffFromCheby = chebyPolyFromZ';

  chebyNumeratorCoeffs = bestChebyCoeffs(1:bestNumTermCount);
  chebyDenominatorCoeffs = [ 1; bestChebyCoeffs(bestNumTermCount+1:end) ];
  zNumeratorCoeffs = zCoeffFromCheby(1:bestNumTermCount,1:bestNumTermCount)*chebyNumeratorCoeffs;
  zDenominatorCoeffs = zCoeffFromCheby(1:bestDenomTermCount, 1:bestDenomTermCount)*chebyDenominatorCoeffs;
  %chebyPolyFromZ
  %chebyNumeratorCoeffs
  %chebyDenominatorCoeffs
  %zNumeratorCoeffs
  %zDenominatorCoeffs
  zNumeratorCoeffs = zNumeratorCoeffs ./ zDenominatorCoeffs(1);
  zDenominatorCoeffs = zDenominatorCoeffs ./ zDenominatorCoeffs(1);
endfunction
  
function [hornNumeratorCoeffs, hornDenominatorCoeffs] = findHornFromZ(zNumeratorCoeffs, zDenominatorCoeffs, lb, ub)
  numCount=length(zNumeratorCoeffs);
  denomCount=length(zDenominatorCoeffs);
  maxTermCount=max(numCount, denomCount);
  
  hornPolyFromZ = zeros(maxTermCount,maxTermCount);
  hornPolyFromZ(1,1) = 1;
  for ii=2:maxTermCount
    hornPolyFromZ(ii,:) = [0 hornPolyFromZ(ii-1, 1:end-1)]*(ub-lb)/2 + hornPolyFromZ(ii-1,:)*(ub+lb)/2;
  endfor

  % combine them into one nice converter.

  if (diag(diag(hornPolyFromZ)) == hornPolyFromZ)
    hornCoeffFromZ = diag(1./diag(hornPolyFromZ));
  else
    oldSetting = warning("query", "singular-matrix-div");
    warning("off", "singular-matrix-div");
    hornCoeffFromZ = inv(hornPolyFromZ');
    warning(oldSetting.state, "singular-matrix-div");
  endif

  hornNumeratorCoeffs = hornCoeffFromZ(1:numCount,1:numCount)*zNumeratorCoeffs;
  hornDenominatorCoeffs = hornCoeffFromZ(1:denomCount, 1:denomCount)*zDenominatorCoeffs;
  %hornPolyFromZ
  %hornCoeffFromZ
  %zNumeratorCoeffs
  %zDenominatorCoeffs
  %hornNumeratorCoeffs
  %hornDenominatorCoeffs
  hornNumeratorCoeffs = hornNumeratorCoeffs ./ hornDenominatorCoeffs(1);
  hornDenominatorCoeffs = hornDenominatorCoeffs ./ hornDenominatorCoeffs(1);
endfunction

function [approximation,coeffs] = buildApproximate(yyy, interpPoly, numTermCount, denomTermCount)

  coeffCount = numTermCount+denomTermCount-1;
  AAA = zeros(rows(yyy), coeffCount);
  bbb = yyy;
  for ii=1:numTermCount
    AAA(:,ii) = interpPoly(:,ii);
  endfor
  for jj=2:denomTermCount
    AAA(:,jj+numTermCount-1) = -yyy.*interpPoly(:,jj);
  endfor

  coeffs = AAA\bbb; % find the coeffs
  approximation = interpPoly(:,1:numTermCount)*coeffs(1:numTermCount) ./ (ones(size(yyy)) + interpPoly(:,2:denomTermCount)*coeffs(numTermCount+1:end));
endfunction
