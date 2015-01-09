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
        
        Jama.Matrix posMatrix = new Matrix(positions, positions.length, positions[0].length);
        Jama.Matrix transposePosMatrix = posMatrix.transpose();
        Jama.Matrix distanceMatrix = new Matrix(distances, distances.length);
        
//        Jama.Matrix posTimesPosTranspose = posMatrix.arrayTimes(transposePosMatrix);
        
        // Calculate a
        for(int i = 0; i < distances.length; i++) {
            Jama.Matrix p_i = new Matrix(positions[i], positions[i].length);
            Jama.Matrix p_i_transpose = p_i.transpose();
            
            Jama.Matrix val = p_i.times(p_i_transpose);
            System.out.println(i);
            System.out.println("p_i");
            p_i.print(3, 3);
            System.out.println("p_i_transpose");
            p_i_transpose.print(3, 3);
            System.out.println("val");
            val.print(3, 3);
        }
//        Jama.Matrix temp1 = posTimesPosTranspose.arrayTimes(posMatrix);
//        Jama.Matrix temp2 = distanceMatrix.arrayTimes(distanceMatrix).arrayTimes(posMatrix);
//        Jama.Matrix a = temp1.minus(temp2);
//        
//        // Calculate B
//        temp1 = posTimesPosTranspose.times(-2);
//        temp2 = posTimesPosTranspose.times(arg0);
        
        return null;
    }

    public static void main(String[] args) {

        TrilaterationModule triMod = new TrilaterationModule();
        
        double[][] positions = {{3., 4., 5.},
                                {6., 7., 8.},
                                {9., 10., 11.}};
        
        double[] distances = {5., 3., 2.};
        
        triMod.calculateLocation(distances, positions);
    }
}