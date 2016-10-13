function [hornNumerator, hornDenominator] = rationalApprox(fff, lb, ub, numPoints, tolerance)
  zzz = linspace(-1,1,numPoints)';
  xxx = zzz*(ub-lb)/2 + (ub+lb)/2;
  yyy = fff(xxx);

  for cost=3:64
    [hornNumerator, hornDenominator, thisTolerance] = findBestAtCost(yyy, zzz, lb, ub, cost);
    numTermCount = rows(hornNumerator);
    denomTermCount = rows(hornDenominator);
    thisTolerance

    if (thisTolerance < tolerance)
      %print it out
      numeratorString = cell(size(hornNumerator));
      denominatorString = cell(size(hornDenominator));
      for ii=1:length(hornNumerator)
        numeratorString{ii} = sprintf("%.17e", hornNumerator(ii));
      endfor
      for ii=1:length(hornDenominator)
        denominatorString{ii} = sprintf("%.17e", hornDenominator(ii));
      endfor
      printf("if (useLookups)");
      printf("{\n");
      printf("   //rationalApprox(%s, %f, %f, %d, %e)\n", func2str(fff), lb, ub, numPoints, tolerance);
      printf("   double x = ;\n");
      printf("   const int pSize = %d;\n", numTermCount);
      printf("   const double pCoeff[] = { %s }\;\n", strjoin(numeratorString, ", "));
      printf("   double numerator=pCoeff[pSize-1];\n");
      printf("   for (int ip=pSize-2; ip>=0; --ip)\n");
      printf("   {\n");
      printf("      numerator = pCoeff[ip] + x*numerator;\n")
      printf("   }\n");
      if (denomTermCount == 1)
        printf("   ret = numerator;\n");
      else
        printf("   const int qSize = %d;\n", denomTermCount);
        printf("   double qCoeff[] = { %s };\n", strjoin(denominatorString, ", "));
        printf("   double denominator = qCoeff[qSize-1];\n");
        printf("   for (int iq=qSize-2; iq>=0; --iq)\n");
        printf("   {\n");
        printf("      denomerator = qCoeff[iq] + x*denominator;\n");
        printf("   }\n");
        printf("   ret = numerator/denominator;\n");
      endif
      printf("}\n");
      printf("else\n");
      printf("{\n\n");
      printf("}\n");
      break;
    endif
  endfor
endfunction

function [hornNumeratorCoeffs, hornDenominatorCoeffs, bestTolerance] = findBestAtCost(yyy, zzz, lb, ub, maxTermCount)
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

  bestChebyCoeffs = zeros(maxTermCount-1);
  bestTolerance = 1e30;
  bestNumTermCount = 0;
  bestDenomTermCount = 0;
  for denomTermCount=1:maxTermCount-1
    numTermCount = maxTermCount - denomTermCount;
    
    [approximation, chebyCoeffs] = buildApproximate(yyy, chebyPoly, numTermCount, denomTermCount);
    tolerance = norm(approximation - yyy)/rows(yyy);
    if (tolerance < bestTolerance)
      bestTolerance = tolerance;
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

  % build the horn coeff
  zFromCheby = zeros(maxTermCount, maxTermCount);
  zFromCheby(1:2,1:2) = eye(2);
  for ii=3:maxTermCount
    zFromCheby(:,ii) = [0; 2*zFromCheby(1:end-1, ii-1)] - zFromCheby(:, ii-2);
  endfor

  zFromHorn = zeros(maxTermCount, maxTermCount);
  zFromHorn(1,1) = 1;
  for ii=2:maxTermCount
    zFromHorn(:,ii) = [0; zFromHorn(1:end-1, ii-1)]*(ub-lb)/2 + zFromHorn(:,ii-1)*(ub+lb)/2;
  endfor

  % combine them into one nice converter.
  oldSetting = warning("query", "singular-matrix-div");
  warning("off", "singular-matrix-div");
  hornFromCheby = zFromHorn \ zFromCheby;
  warning(oldSetting.state, "singular-matrix-div");
  
  chebyNumeratorCoeffs = bestChebyCoeffs(1:bestNumTermCount);
  chebyDenominatorCoeffs = [ 1; bestChebyCoeffs(bestNumTermCount+1:end) ];
  hornNumeratorCoeffs = hornFromCheby(1:bestNumTermCount,1:bestNumTermCount)*chebyNumeratorCoeffs;
  hornDenominatorCoeffs = hornFromCheby(1:bestDenomTermCount, 1:bestDenomTermCount)*chebyDenominatorCoeffs;
  %zFromCheby
  %zFromHorn
  %hornFromCheby
  %chebyNumeratorCoeffs
  %chebyDenominatorCoeffs
  %hornNumeratorCoeffs
  %hornDenominatorCoeffs
  hornNumeratorCoeffs = hornNumeratorCoeffs ./ hornDenominatorCoeffs(1);
  hornDenominatorCoeffs = hornDenominatorCoeffs ./ hornDenominatorCoeffs(1);
  hornCoeffs = [hornNumeratorCoeffs; hornDenominatorCoeffs(2:end)];

  
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
  approximation = AAA*coeffs;
endfunction
