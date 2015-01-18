import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.linear.QRDecomposition;
import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.apache.commons.math3.linear.MatrixUtils;

public class TrilaterationModule {

    /**
     * Class constructor.
     */
    public TrilaterationModule() {

    }

    public RectCoordinates calculateLocation(double[] distances, double[][] positions) {
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
        
        if(numDimensions == 3) {
            
        }
        else if(numDimensions == 2) {
            
        }
                
        System.out.println("a");
        printRealMatrix(a);
        System.out.println("b");
        printRealMatrix(b);
        System.out.println("c");
        printRealMatrix(c);
        System.out.println("f");
        printRealMatrix(f);
        System.out.println("fPrime");
        printRealMatrix(fPrime);
        System.out.println("h");
        printRealMatrix(h);
        System.out.println("hPrime");
        printRealMatrix(hPrime);
        System.out.println("q");
        printRealMatrix(q);
        System.out.println("u");
        printRealMatrix(u);
        System.out.println("qT*q = " + qTransposeTimesQ);
        
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

    public static void main(String[] args) {

        TrilaterationModule triMod = new TrilaterationModule();
        
        double[][] positions = {{0., 0., 0.},
                                {0., 2., 2.},
                                {3., 0., 4.},
                                {2., 5., 1.}};
        
        double[] distances = {4.2426, 3.7417, 5., 2.4495};
        
        triMod.calculateLocation(distances, positions);
    }
}