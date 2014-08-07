{
    p.boundaryField().updateCoeffs();

    volScalarField rAU(1.0/UEqn().A());
    U = rAU*UEqn().H();
    UEqn.clear();

    phi = fvc::interpolate(U, "interpolate(HbyA)") & mesh.Sf();
    adjustPhi(phi, U, p);

    // Non-orthogonal pressure corrector loop
    for (int nonOrth=0; nonOrth<=simple.nNonOrthCorr(); nonOrth++)
    {
        fvScalarMatrix pEqn
        (
            fvm::laplacian(rAU, p) == fvc::div(phi)
        );

        pEqn.setReference(pRefCell, pRefValue);

        pEqn.solve();
//-----------------------------------------------               
               List<List<scalar> > A; 
               List<scalar> b;
        Info << "Matrixes created!" << nl << endl;        
                double upSum1 = 0;			//
                double downSum1 = 0;			//
                double upSum2 = 0;
		double downSum2 = 0;
                double F1 = 0;				//
		double F2 = 0;
		
                A.resize(p.size());			// A and p have equal size
                b.resize(p.size());
Info << "Matrixes resized!" << nl << endl; 
		forAll(A, i)
                {
            	    A[i].resize(p.size());		// every element of A is array
            	    forAll(A[i],j)			// clearing A and b
            	    {
            		A[i][j] = 0.0;
            	    }
            	    b[i] = 0.0;
                }
Info << "Matrixes cleared!" << nl << endl; 
// //                 arma::mat Ar = arma::zeros<arma::mat>(p.size(),p.size());
// // 		arma::mat br = arma::zeros<arma::mat>(1,p.size());

//                 forAll(p,i) // forAll(A,i)
//                 {
//             	    Ar(i,i) = pEqn.diag()[i];
// 		    br(i) = pEqn.source()[i];
// 		}
		forAll(p,i)
		{
		  A[i][i] = pEqn.diag()[i];
		  b[i] = pEqn.source()[i];
		}
Info << "Diag and b filled!" << nl << endl; 
                const lduAddressing& addr = pEqn.lduAddr();
                const labelList& lowerAddr = addr.lowerAddr();
                const labelList& upperAddr = addr.upperAddr();                
                
                forAll(lowerAddr, i)
                {
// 		    Ar(lowerAddr[i],upperAddr[i]) = pEqn.upper()[i];
// 		    Ar(upperAddr[i],lowerAddr[i]) = pEqn.lower()[i];  
		    A[lowerAddr[i],upperAddr[i]] = pEqn.upper()[i];
		    A[upperAddr[i],lowerAddr[i]] = pEqn.lower()[i];
                }
Info << "Upper Lower filled!" << nl << endl; 

                forAll(p.boundaryField(),I)
                {
            	    const fvPatch &ptch=p.boundaryField()[I].patch();
            	    forAll(ptch,J)
            	    {
           		int w=ptch.faceCells()[J];
// 			Ar(w,w)+=pEqn.internalCoeffs()[I][J];
// 			br(w) +=pEqn.boundaryCoeffs()[I][J];
			A[w][w]+=pEqn.internalCoeffs()[I][J];
			b[w] +=pEqn.boundaryCoeffs()[I][J];
            	    }
                
                }
                
Info << "Boundary filled!" << nl << endl; 

/*                arma::mat Ar2 = Ar;
		arma::inplace_trans(Ar);
		arma::mat G = Ar*Ar2;

                arma::mat D = arma::zeros<arma::mat>(p.size(),p.size()); */// = pow(Ar%Eye,-0.5); //CHIT!!

// 		for (int i=0; i<p.size(); i++) // forAll(A,i)
//                 {
// 		  for (int j=0; j<p.size(); j++)
// 		  {
// 		    if (i==j)
// 		    {
// 			D(i,j) = 1/std::pow(G(i,j),0.5);
// 		    }
// 		  }
// 		}
		
// 		arma::mat R = D*G*D;
		upSum1 = 0;
		upSum2 = p.size();
		
// 		for (int i=0; i<p.size(); i++) // forAll(A,i)
//                 {
// 		  for (int j=0; j<p.size(); j++)
// 		  {
// 		    if (i==j)
// 		    {
// 		      upSum1 += R(i,j); ////
// 		    }
// 
// 		    downSum1 += R(i,j)*R(i,j);
// 		  }
// 		}
Info << "Out!" << nl << endl; 
		Info << "getR=" << getR(1,1,p.size(), A) << nl << endl;
		Info << "getR=" << getG(1,1,p.size(), A) << nl << endl;
		Info << "getR=" << getD(1,1,p.size(), A) << nl << endl;
				
                F1 = upSum1*upSum1 / downSum1;
                F2 = upSum2*upSum2 / downSum2;
		
		Info << "upSum1=" << upSum1 << nl << endl;
		Info << "downSum1=" << downSum1 << nl << endl;
		Info << "F1=" << F1 << nl << endl;
		Info << "upSum2=" << upSum2 << nl << endl;
		Info << "downSum2=" << downSum2 << nl << endl;
		Info << "F2=" << F2 << nl << endl;
		

		
//----------------------------------------------------------------- 
        if (nonOrth == simple.nNonOrthCorr())
        {
            phi -= pEqn.flux();
        }
    }

    #include "continuityErrs.H"

    // Explicitly relax pressure for momentum corrector
    p.relax();

    // Momentum corrector
    U -= rAU*fvc::grad(p);
    U.correctBoundaryConditions();

  
}

