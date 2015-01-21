package walker.blue.tri.lib;

import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.linear.QRDecomposition;
import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.apache.commons.math3.linear.MatrixUtils;

public class TrilaterationModule {
    private RealMatrix previousPosition;
    /**
     * Class constructor.
     */
    public TrilaterationModule() {
        
    }

    public double[][] calculateLocation(double[] distances, double[][] positions) {
        int numDimensions = positions[0].length;
        int numSamples = positions.length;
         
        RealMatrix a = new Array2DRowRealMatrix(numDimensions, 1);
        RealMatrix b = new Array2DRowRealMatrix(numDimensions, numDimensions);
        RealMatrix c = new Array2DRowRealMatrix(numDimensions, 1);
        RealMatrix h = new Array2DRowRealMatrix(numDimensions, numDimensions);
        double qTransposeTimesQFirstTerm = 0;
        double qTransposeTimesQSecondTerm = 0;
        // Calculate a, b, c, and partially calculate h
        for(int i = 0; i < numSamples; i++) {
            RealMatrix position = new Array2DRowRealMatrix(positions[i]);
            RealMatrix transposePosition = position.transpose();
            RealMatrix posTimesTranspose = position.multiply(transposePosition);
            double distanceSquared = distances[i]*distances[i];
            qTransposeTimesQSecondTerm += distanceSquared;
            
            // a
            RealMatrix aLeftTerm = posTimesTranspose.multiply(position);
            RealMatrix aRightTerm = position.scalarMultiply(distanceSquared);
            a = a.add(aLeftTerm.subtract(aRightTerm));
  
            // b
            RealMatrix bLeftTerm = posTimesTranspose.scalarMultiply(-2);
            RealMatrix identity = MatrixUtils.createRealIdentityMatrix(numDimensions);
            double transposeTimesPos = transposePosition.multiply(position).getEntry(0, 0);
            RealMatrix bMiddleTerm = identity.scalarMultiply(transposeTimesPos);
            RealMatrix bRightTerm = identity.scalarMultiply(distanceSquared);
            b = b.add(bLeftTerm.subtract(bMiddleTerm).add(bRightTerm));
                    
            // c
            c = c.add(position);
            
            // h
            h = h.add(posTimesTranspose);
            
            // qTransposeTimesQFirstTerm
            qTransposeTimesQFirstTerm += transposeTimesPos;
        }
        
        double inverseNumSamples = 1.0/numSamples;
        a = a.scalarMultiply(inverseNumSamples); // a's calculation done
        b = b.scalarMultiply(inverseNumSamples); // b's calculation done
        c = c.scalarMultiply(inverseNumSamples); // c's calculation done
        
        RealMatrix twoCTimesCTranspose = c.scalarMultiply(2).multiply(c.transpose());
        // Calculate f
        RealMatrix f = a.add(b.multiply(c)).add(twoCTimesCTranspose.multiply(c)); // f's calculation done
        
        // Calculate fPrime
        RealMatrix fPrime = new Array2DRowRealMatrix(numDimensions - 1, 1);
        for(int i = 0; i < fPrime.getRowDimension(); i++) {
            fPrime.setEntry(i, 0, f.getEntry(i, 0) - f.getEntry(numDimensions - 1, 0));
        }
        
        h = h.scalarMultiply(-2.0/numSamples).add(twoCTimesCTranspose); // h's calculation done
        
        // Calculate hPrime
        RealMatrix hPrime = new Array2DRowRealMatrix(numDimensions - 1, numDimensions);
        for(int i = 0; i < hPrime.getRowDimension(); i++) {
            for(int j = 0; j < hPrime.getColumnDimension(); j++) {
                hPrime.setEntry(i, j, h.getEntry(i, j) - h.getEntry(numDimensions - 1, j));
            }
        }
        
        QRDecomposition qrDecomp = new QRDecomposition(hPrime);
        RealMatrix q = qrDecomp.getQ();
        RealMatrix u = qrDecomp.getR();
        
        qTransposeTimesQFirstTerm *= -inverseNumSamples;
        qTransposeTimesQSecondTerm /= numSamples;
        double qTransposeTimesQThirdTerm = c.transpose().multiply(c).getEntry(0, 0);
        double qTransposeTimesQ = qTransposeTimesQFirstTerm + qTransposeTimesQSecondTerm + qTransposeTimesQThirdTerm;
        
        RealMatrix v = q.transpose().multiply(fPrime);
        
        if(numDimensions == 3) {
            double u11 = u.getEntry(0, 0);
            double u12 = u.getEntry(0, 1);
            double u13 = u.getEntry(0, 2);
            double u22 = u.getEntry(1, 1);
            double u23 = u.getEntry(1, 2);
            double v1 = v.getEntry(0, 0);
            double v2 = v.getEntry(1, 0);
            // Solving for values a, b, and c to be passed into the quadratic formula
            double quadraticA = Math.pow(((u12*u23)/(u11*u22)) - (u13/u11), 2) + Math.pow(u23/u22, 2) + 1;         
            double quadraticB = 2 * ((((u12*v2)/(u11*u22))-( v1/u11 )) * (((u12*u23)/(u11*u22))-(u13/u11)) + ((u23*v2)/(u22*u22)));                       
            double quadraticC = Math.pow(((u12*v2)/(u11*u22))-(v1/u11), 2) + Math.pow(v2/u22, 2) - qTransposeTimesQ;
            
            // Using the quadratic equation, two possible solutions for q3 can be found
            double sqrtTerm = Math.sqrt((Math.pow(quadraticB, 2)) - (4*quadraticA*quadraticC));
            double twoA = 2*quadraticA;
            
            RealMatrix qSol1 = new Array2DRowRealMatrix(numDimensions, 1);
            RealMatrix qSol2 = new Array2DRowRealMatrix(numDimensions, 1);
            
            // First q solution
            qSol1.setEntry(2, 0, (-quadraticB + sqrtTerm)/twoA);
            qSol1.setEntry(0, 0, threeDimensionalSolveQ1(u, v, qSol1.getEntry(2, 0)));
            qSol1.setEntry(1, 0, threeDimensionalSolveQ2(u, v, qSol1.getEntry(2, 0)));
            // Second q solution
            qSol2.setEntry(2, 0, (-quadraticB - sqrtTerm)/twoA);
            qSol2.setEntry(0, 0, threeDimensionalSolveQ1(u, v, qSol2.getEntry(2, 0)));
            qSol2.setEntry(1, 0, threeDimensionalSolveQ2(u, v, qSol2.getEntry(2, 0)));
            
            RealMatrix pSol1 = qSol1.add(c);
            RealMatrix pSol2 = qSol2.add(c);
            
            System.out.println("pSol1");
            printRealMatrix(pSol1);
            System.out.println("pSol2");
            printRealMatrix(pSol2);
            
            RealMatrix solution;
            
            // If this is the first calculated position, then go through
            // the longer process of using the distance formula to determine
            // which of the solutions is most accurate. Otherwise, if it
            // is not the first calculated position, then just return the 
            // solution that is closest to the previous position.
            if(previousPosition == null) {
                System.out.println("First calculation");
                solution = calculateBestSolution(distances, positions, pSol1, pSol2);
            }
            else {
                System.out.println("After first calculation");
                solution = solutionClosestToPrevious(pSol1, pSol2);
            }
            
            previousPosition = solution;
            
            return solution.getData();
        }
        else if(numDimensions == 2) {
            // TODO: Write calculations for two-dimensions
        }        
        return null;
    }
    
    public void printRealMatrix(RealMatrix matrix) {
        for(int i = 0; i < matrix.getRowDimension(); i++) {
            for(int j = 0; j < matrix.getColumnDimension(); j++) {
                System.out.print(matrix.getEntry(i, j) + " ");
            }
            System.out.println();
        }
        System.out.println();
    }
    
    private double threeDimensionalSolveQ1(RealMatrix u, RealMatrix v, double q3) {
        double u11 = u.getEntry(0, 0);
        double u12 = u.getEntry(0, 1);
        double u13 = u.getEntry(0, 2);
        double u22 = u.getEntry(1, 1);
        double u23 = u.getEntry(1, 2);
        double v1 = v.getEntry(0, 0);
        double v2 = v.getEntry(1, 0);
        
        return (((u12*v2)/(u11*u22)) - (v1/u11)) + ((((u12*u23)/(u11*u22)) - (u13/u11))*q3);
    }
    
    private double threeDimensionalSolveQ2(RealMatrix u, RealMatrix v, double q3) {
        double u22 = u.getEntry(1, 1);
        double u23 = u.getEntry(1, 2);
        double v2 = v.getEntry(1, 0);
        
        return (-v2/u22) - ((u23/u22)*q3);
    }
    
    private RealMatrix calculateBestSolution(double[] distances, double[][] positions, RealMatrix solutionA, RealMatrix solutionB) {
        // Increments whenever the particular solutions position is more accurate then the
        // other for a given inputed position.
        int countA = 0;
        int countB = 0;
        // Compare solutionA's position to solutionB's position by calculating the
        // distances between each solution's position and each of the inputed positions
        // to see which one is closest to the inputed distance measurements
        for(int i = 0; i < positions.length; i++) {
            double solADistance = Math.sqrt(Math.pow(solutionA.getEntry(0, 0) - positions[i][0], 2) +
                    Math.pow(solutionA.getEntry(1, 0) - positions[i][1], 2) +
                    Math.pow(solutionA.getEntry(2, 0) - positions[i][2], 2));
            double solBDistance = Math.sqrt(Math.pow(solutionB.getEntry(0, 0) - positions[i][0], 2) +
                    Math.pow(solutionB.getEntry(1, 0) - positions[i][1], 2) +
                    Math.pow(solutionB.getEntry(2, 0) - positions[i][2], 2));
            
            double solADiff = solADistance - distances[i];
            double solBDiff = solBDistance - distances[i];
            // To prevent negative differences
            if(solADiff < 0) {
                solADiff = -solADiff;
            }
            if(solBDiff < 0) {
                solBDiff = -solBDiff;
            }
            
            if(solADiff <= solBDiff) {
                countA++;
            }
            else {
                countB++;
            }
        }
        
        if(countA >= countB) {
            return solutionA;
        }
        else {
            return solutionB;
        }
    }
    
    private RealMatrix solutionClosestToPrevious(RealMatrix solutionA, RealMatrix solutionB) {
        double solADistance = Math.sqrt(Math.pow(solutionA.getEntry(0, 0) - previousPosition.getEntry(0, 0), 2) +
                Math.pow(solutionA.getEntry(1, 0) - previousPosition.getEntry(1, 0), 2) +
                Math.pow(solutionA.getEntry(2, 0) - previousPosition.getEntry(2, 0), 2));
        double solBDistance = Math.sqrt(Math.pow(solutionB.getEntry(0, 0) - previousPosition.getEntry(0, 0), 2) +
                Math.pow(solutionB.getEntry(1, 0) - previousPosition.getEntry(1, 0), 2) +
                Math.pow(solutionB.getEntry(2, 0) - previousPosition.getEntry(2, 0), 2));
        
        if(solADistance <= solBDistance) {
            return solutionA;
        }
        else {
            return solutionB;
        }
    }

    public static void main(String[] args) {

        TrilaterationModule triMod = new TrilaterationModule();
        
        double[][] positions = {{0., 0., 0.},
                                {0., 2., 2.},
                                {3., 0., 4.},
                                {2., 5., 1.}};
        
        double[] distances1 = {4.1231056256177, 3, 2.8284271247462, 3.3166247903554};
        double[] distances2 = {3., 3., 4.5, 4};
        double[] distances3 = {4.2426, 3.7417, 5, 2.4495};
        double[] distances4 = {38.858, 36.689, 36.125, 33.615};
        double[][] calcedLocation = triMod.calculateLocation(distances4, positions);
        for(int i = 0; i < calcedLocation.length; i++) {
            System.out.println(calcedLocation[i][0]);
        }
        calcedLocation = triMod.calculateLocation(distances1, positions);
        for(int i = 0; i < calcedLocation.length; i++) {
            System.out.println(calcedLocation[i][0]);
        }
        calcedLocation = triMod.calculateLocation(distances1, positions);
        for(int i = 0; i < calcedLocation.length; i++) {
            System.out.println(calcedLocation[i][0]);
        }
    }
}
