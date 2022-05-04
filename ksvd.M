function [ mD ] = DL_KSVD( mD, mX, p )

      numberOfAtoms = size(mD, 2);

      for ii = 1:50
          %<! Update the Representations
          mA = algo_omp(mX, mD, p);

          %<! Update the Dictionary
          for jj = 1:5
              mE = mX - (mD * mA);
              for kk = 1:numberOfAtoms
                  vP         = find(mA(kk, :));
                  mEP        = mE(:, vP) + mD(:, kk) * mA(kk, vP);
                  vA         = mD(:, kk) * mEP;
                  mA(kk, vP) = vA;
                  mD(:, kk)  = mEP * vA' / (vA * vA');
              end
          end

          mD = bsxfun(@rdivide, mD, sqrt(sum(mD .^ 2, 1)));
          
      end

end
