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
        
        Jama.Matrix a = new Matrix(positions[0].length, 1);
        Jama.Matrix b = new Matrix(positions.length, positions[0].length);
        Jama.Matrix c = new Matrix(positions[0].length, 1);
        // Calculate a and B and c
        for(int i = 0; i < positions.length; i++) {
            Jama.Matrix position = new Matrix(positions[i], positions[i].length);
            Jama.Matrix transposePosition = position.transpose();
            Jama.Matrix posTimesTranspose = position.times(transposePosition);
            double distanceSquared = distances[i]*distances[i];
            
            // a
            Jama.Matrix aLeftTerm = posTimesTranspose.times(position);
            Jama.Matrix aRightTerm = position.times(distanceSquared); 
            a.plusEquals(aLeftTerm.minus(aRightTerm));
  
            // B
            Jama.Matrix bLeftTerm = posTimesTranspose.times(-2);
            //Jama.Matrix bMiddleTerm = ;
            //Jama.Matrix bRightTerm = ;
                    
            // c
            c.plusEquals(position);    
        }
        
        a.times(1.0/positions.length); // a's calculation done
        b.times(1.0/positions.length); // b's calculation done
        c.times(1.0/positions.length); // c's calculation done
        
        System.out.println("a");
        a.print(3, 3);
        
        System.out.println("c");
        c.print(3, 3);
        
        return null;
    }

    public static void main(String[] args) {

        TrilaterationModule triMod = new TrilaterationModule();
        
        double[][] positions = {{0., 1., 2.},
                                {3., 4., 5.},
                                {6., 7., 8.}};
        
        double[] distances = {5., 3., 2.};
        
        triMod.calculateLocation(distances, positions);
    }
}