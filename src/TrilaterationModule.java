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
        // Calculate a     
        for(int i = 0; i < distances.length; i++) {
            Jama.Matrix position = new Matrix(positions[i], positions[i].length);
            Jama.Matrix transposePosition = position.transpose();
            
            Jama.Matrix leftTerm = position.times(transposePosition).times(position);
            Jama.Matrix rightTerm = position.times(distances[i]*distances[i]);
            
            a.plusEquals(leftTerm.minus(rightTerm));
        }
        
        System.out.println("a");
        a.print(3, 3);
        
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