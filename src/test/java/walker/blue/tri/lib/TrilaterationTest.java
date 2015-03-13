package walker.blue.tri.lib;

import org.junit.Assert;
import org.junit.Before;
import org.junit.Test;

public class TrilaterationTest {

    Trilateration trilateration;

    public static double error = 0.001;

    @Before
    public void classSetup() {
        trilateration = new Trilateration();
    }

    @Test
    public void test1() {
        double[][] positions = {{0., 0., 0.},
                {0., 2., 2.},
                {3., 0., 4.},
                {2., 5., 1.}};

        double[][] exactLocation = {{-12.0},{-1.0},{-2.0}};

        double[] distances = calculateDistances(exactLocation, positions);

        double[][] calculatedLocation = trilateration.calculateLocation(distances, positions);

        Assert.assertTrue(isWithinErrorRadius(exactLocation, calculatedLocation));
    }

    @Test
    public void test2() {
        double[][] positions = {{0., 0., 0.},
                {0., 2., 2.},
                {3., 0., 4.},
                {2., 5., 1.}};

        double[][] exactLocation = {{1002.101},{10.0},{-223.00988889}};

        double[] distances = calculateDistances(exactLocation, positions);

        double[][] calculatedLocation = trilateration.calculateLocation(distances, positions);

        Assert.assertTrue(isWithinErrorRadius(exactLocation, calculatedLocation));
    }

    @Test
    public void test3() {
        double[][] positions = {{0., 0., 0.},
                {0., 2., 2.},
                {3., 0., 4.},
                {2., 5., 1.}};

        double[][] exactLocation = {{1.111},{0.000001},{3.456}};

        double[] distances = calculateDistances(exactLocation, positions);

        double[][] calculatedLocation = trilateration.calculateLocation(distances, positions);

        Assert.assertTrue(isWithinErrorRadius(exactLocation, calculatedLocation));
    }

    @Test
    public void test4() {
        double[][] positions = {{0., 0., 0.},
                {0., 2., 2.},
                {3., 0., 4.},
                {2., 5., 1.}};

        double[][] exactLocation = {{24.333},{10.0000001},{-2.0}};

        double[] distances = calculateDistances(exactLocation, positions);

        double[][] calculatedLocation = trilateration.calculateLocation(distances, positions);

        Assert.assertTrue(isWithinErrorRadius(exactLocation, calculatedLocation));
    }

    private double[] calculateDistances(double[][] exactLocation, double[][] otherLocations) {
        double[] distances = new double[otherLocations.length];

        for(int i = 0; i < distances.length; i++) {
            distances[i] = Math.sqrt(
                    Math.pow(exactLocation[0][0] - otherLocations[i][0], 2) +
                            Math.pow(exactLocation[1][0] - otherLocations[i][1], 2) +
                            Math.pow(exactLocation[2][0] - otherLocations[i][2], 2)
            );
        }

        return distances;
    }

    private boolean isWithinErrorRadius(double[][] exactLocation, double[][] calculatedLocation) {
        // If the error of the calculatedLocation in one direction is
        // greater than the allowed error, then return false.
        // Otherwise, return true.
        for(int i = 0; i < exactLocation.length; i++) {
            double errorAlongAxis = calculateError(exactLocation[i][0], calculatedLocation[i][0]);
            if(errorAlongAxis > this.error) {
                return false;
            }
        }
        return true;
    }

    private double calculateError(double exactVal, double calculatedVal) {
        return Math.abs(calculatedVal - exactVal) / Math.abs(exactVal);
    }
}
