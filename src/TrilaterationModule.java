import java.util.List;

import Jama.Matrix;
import Jama.QRDecomposition;

public class TrilaterationModule {

    /**
     * Class constructor.
     */
    public TrilaterationModule() {

    }

    public RectCoordinates calculateLocation(double[] distances, double[][] positions) {
        int numDimensions = positions[0].length;
        int numSamples = positions.length;
        
        Jama.Matrix a = new Matrix(numDimensions, 1);  
        Jama.Matrix b = new Matrix(numDimensions, numDimensions); 
        Jama.Matrix c = new Matrix(numDimensions, 1);
        Jama.Matrix h = new Matrix(numDimensions, numDimensions);
        // Calculate a, b, c, and partially calculate h
        for(int i = 0; i < numSamples; i++) {
            Jama.Matrix position = new Matrix(positions[i], numDimensions);
            Jama.Matrix transposePosition = position.transpose();
            Jama.Matrix posTimesTranspose = position.times(transposePosition);
            double distanceSquared = distances[i]*distances[i];
            
            // a
            Jama.Matrix aLeftTerm = posTimesTranspose.times(position);
            Jama.Matrix aRightTerm = position.times(distanceSquared); 
            a.plusEquals(aLeftTerm.minus(aRightTerm));
  
            // b
            Jama.Matrix bLeftTerm = posTimesTranspose.times(-2);
            Jama.Matrix identity = Jama.Matrix.identity(numDimensions, numDimensions);
            Jama.Matrix bMiddleTerm = identity.times(transposePosition.times(position).get(0, 0));
            Jama.Matrix bRightTerm = identity.times(distanceSquared);
            b.plusEquals(bLeftTerm.minus(bMiddleTerm).plus(bRightTerm));
                    
            // c
            c.plusEquals(position);
            
            // h
            h.plusEquals(posTimesTranspose);
        }
        
        double inverseNumSamples = 1.0/numSamples;
        a.times(inverseNumSamples); // a's calculation done
        b.times(inverseNumSamples); // b's calculation done
        c.times(inverseNumSamples); // c's calculation done
         
        Jama.Matrix twoCTimesCTranspose = c.times(2).times(c.transpose());
        // Calculate f
        Jama.Matrix f = a.plus(b.times(c)).plus(twoCTimesCTranspose.times(c)); // f's calculation done
        
        // Calculate fPrime
        Jama.Matrix fPrime = new Matrix(numDimensions - 1, 1);
        for(int i = 0; i < fPrime.getRowDimension(); i++) {
            fPrime.set(i, 0, f.get(i, 0) - f.get(numDimensions - 1, 0));
        }
        
        h.times(-2.0/numSamples).plus(twoCTimesCTranspose); // f's calculation done
        
        // Calculate hPrime
        Jama.Matrix hPrime = new Matrix(numDimensions - 1, numDimensions);
        for(int i = 0; i < hPrime.getRowDimension(); i++) {
            for(int j = 0; j < hPrime.getColumnDimension(); j++) {
                hPrime.set(i, j, h.get(i, j) - h.get(numDimensions - 1, j));
            }
        }
        
        //QRDecomposition hPrimeQRDecomposition = new QRDecomposition(hPrime);
        //Jama.Matrix q = hPrimeQRDecomposition.getQ();
        //Jama.Matrix u = hPrimeQRDecomposition.getR();        
              
        System.out.println("a");
        a.print(3, 3);
        System.out.println("b");
        b.print(3, 3);
        System.out.println("c");
        c.print(3, 3);
        System.out.println("f");
        f.print(3, 3);
        System.out.println("fPrime");
        fPrime.print(3, 3);
        System.out.println("h");
        h.print(3, 3);
        System.out.println("hPrime");
        hPrime.print(3, 3);
        //System.out.println("q");
        //q.print(3, 3);
        //System.out.println("u");
        //u.print(3, 3);
        
        return null;
    }

    public static void main(String[] args) {

        TrilaterationModule triMod = new TrilaterationModule();
        
        double[][] positions = {{0., 1., 2.},
                                {3., 4., 5.},
                                {6., 7., 8.},
                                {9., 10., 11.}};
        
        double[] distances = {1., 1., 1., 1.};
        
        triMod.calculateLocation(distances, positions);
    }
}