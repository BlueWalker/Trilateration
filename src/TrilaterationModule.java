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

//        Jama.Matrix matrix = new Matrix(values);
//        QRDecomposition qrDecomposition = new QRDecomposition(matrix);
//        Jama.Matrix orthogonal = qrDecomposition.getQ();
        int numDimensions = positions[0].length;
        int numSamples = positions.length;
        
        Jama.Matrix a = new Matrix(numDimensions, 1);  
        Jama.Matrix b = new Matrix(numDimensions, numDimensions); 
        Jama.Matrix c = new Matrix(numDimensions, 1);
        // Calculate a, b, and c
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
        }
        
        double inverseNumSamples = 1.0/numSamples;
        a.times(inverseNumSamples); // a's calculation done
        b.times(inverseNumSamples); // b's calculation done
        c.times(inverseNumSamples); // c's calculation done
           
        // Calculate f
        Jama.Matrix f = a.plus(b.times(c)).plus(c.times(2).times(c.transpose()).times(c)); // f's calculation done
        
        // Calculate fPrime
        Jama.Matrix fPrime = new Matrix(numDimensions - 1, 1);
        for(int i = 0; i < fPrime.getRowDimension(); i++) {
            fPrime.set(i, 0, f.get(i, 0) - f.get(numDimensions - 1, 0));
        }
        
              
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