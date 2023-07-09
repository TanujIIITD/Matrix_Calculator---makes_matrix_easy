
import java.util.Random;
import java.util.*;

interface common_methods {
    int determinant(double[][] mat_1);

    double[][] transpose(double[][] mat_1);

}

class matrix implements common_methods {
    public double[][] matrix_generator(int rows, int columns) {
        Scanner scn = new Scanner(System.in);

        double[][] matrix = new double[rows][columns];

        for (int i = 0; i < rows; i++) {
            for (int j = 0; j < columns; j++) {
                System.out.print("enter element " + (i + 1) + "x" + (j + 1) + " : ");
                int ele = scn.nextInt();
                matrix[i][j] = ele;
            }
        }
        System.out.println();
        matrix_print(matrix);
        return matrix;
    }

    public void matrix_print(double[][] matrix) {
        int rows = matrix.length;
        int columns = matrix[0].length;
        for (int i = 0; i < rows; i++) {
            for (int j = 0; j < columns; j++) {
                System.out.print(matrix[i][j] + " ");
            }
            System.out.println();
        }
        System.out.println();
    }

    public double[][] addition(double[][] mat_1, double[][] mat_2) {
        if ((mat_1.length != mat_2.length) || (mat_1[0].length != mat_2[0].length)) {
            System.out.println("error: matrices of different dimensions can not be added");
            double[][] result = new double[0][0];
            result[0][0] = 0;
            return result;
        } else {
            int rows = mat_1.length;
            int columns = mat_1[0].length;
            double[][] result = new double[rows][columns];
            for (int i = 0; i < rows; i++) {
                for (int j = 0; j < columns; j++) {
                    result[i][j] = mat_1[i][j] + mat_2[i][j];
                }
            }
            matrix_print(result);
            return result;
        }
    }

    public double[][] subtraction(double[][] mat_1, double[][] mat_2) {
        if ((mat_1.length != mat_2.length) || (mat_1[0].length != mat_2[0].length)) {
            System.out.println("error: matrices of different dimensions can not be subtracted");
            double[][] result = new double[0][0];
            result[0][0] = 0;
            return result;
        } else {
            int rows = mat_1.length;
            int columns = mat_1[0].length;
            double[][] result = new double[rows][columns];
            for (int i = 0; i < rows; i++) {
                for (int j = 0; j < columns; j++) {
                    result[i][j] = mat_1[i][j] - mat_2[i][j];
                }
            }
            matrix_print(result);
            return result;
        }
    }

    public double[][] multiplication(double[][] mat_1, double[][] mat_2) {
        if (mat_1[0].length != mat_2.length) {
            System.out.println("error: wrong matrix dimensions for multiplication");
            double[][] result = new double[0][0];
            result[0][0] = 0;
            return result;
        }

        double[][] result = new double[mat_1.length][mat_2[0].length];
        for (int i = 0; i < mat_1.length; i++) {
            for (int j = 0; j < mat_2[0].length; j++) {
                for (int f = 0; f < mat_1[0].length; f++) {
                    result[i][j] += mat_1[i][f] * mat_2[f][j];
                }
            }
        }
        matrix_print(result);
        return result;

    }

    public double[][] division(double[][] mat_1, double[][] mat_2) {
        int rows = mat_2.length;
        int columns = mat_2[0].length;

        double[][] mat_2_inverse = new double[rows][columns];
        mat_2_inverse = inverse(mat_2);
        double[][] result = new double[rows][columns];
        result = multiplication(mat_1, mat_2_inverse);
        return result;
    }

    public double[][] element_wise_multiplication(double[][] mat_1, double[][] mat_2) {
        int rows = mat_1.length;
        int columns = mat_1[0].length;
        double[][] result = new double[rows][columns];
        for (int i = 0; i < rows; i++) {
            for (int j = 0; j < columns; j++) {
                result[i][j] = mat_1[i][j] * mat_2[i][j];
            }
        }
        matrix_print(result);
        return result;
    }

    public void element_wise_division(double[][] mat_1, double[][] mat_2) {
        int rows = mat_1.length;
        int columns = mat_1[0].length;
        double[][] result = new double[rows][columns];
        for (int i = 0; i < rows; i++) {
            for (int j = 0; j < columns; j++) {
                if ((int) mat_2[i][j] == 0) {
                    System.out.println("error: an element can not be divisible by zero.");
                    return;
                }
                result[i][j] = mat_1[i][j] / mat_2[i][j];
            }
        }

        matrix_print(result);
    }

    public double[] eigen_values(double[][] asd) {
        int rows = asd.length;
        double[] arr_eigen = new double[rows];

        if (rows == 1) {
            System.out.println("Eigen values is: ");
            System.out.println(asd[0][0]);
            return arr_eigen;
        } else if (rows == 2) {
            double cal = Math.sqrt((traces(asd) * traces(asd)) - 4 * determinant(asd));
            arr_eigen[0] = (traces(asd) + cal) / 2;
            arr_eigen[1] = (traces(asd) - cal) / 2;

            System.out.println("Eigen values are: ");
            System.out.println(arr_eigen[0] + " and " + arr_eigen[1]);
            return arr_eigen;
        } else if (rows == 3) {
            System.out.println("Equation will be");
            System.out.println("(-x3) + (" + traces(asd) + ")(-x2) + ("
                    + (cofactor_help(asd)[0][0] + cofactor_help(asd)[1][1] + cofactor_help(asd)[2][2]) + ")(x) + ("
                    + determinant(asd) + ")");
        }
        return arr_eigen;

    }

    public double[][] cofactor_help(double[][] asd) {
        double[][] cof = new double[3][3];
        for (int i = 0; i < 3; i++) {
            for (int j = 0; j < 3; j++) {
                cof[i][j] = cofactor(asd, i, j);
            }
        }
        return cof;
    }

    public double cofactor(double[][] asd, int row, int col) {
        int rows = asd.length;
        int columns = asd[0].length;
        double[][] mat = new double[rows][columns];
        int r = 0;
        int c = 0;
        for (int i = 0; i < rows; i++) {
            if (i == row) {
                continue;
            }
            for (int j = 0; j < columns; j++) {
                if (j == col) {
                    continue;
                }
                mat[r][c] = asd[i][j];
                c++;
            }
            c = 0;
            r++;
        }
        return Math.pow(-1, row + col) * determinant(mat);
    }

    public double traces(double[][] asd) {
        int rows = asd.length;
        double use = 0;
        for (int i = 0; i < rows; i++) {
            use += asd[i][i];
        }
        return use;
    }

    public int determinant(double[][] mat_1) {
        int result = 0;
        if (mat_1.length != mat_1[0].length) {
            System.out.println("error: matrix must be square matrix for determinant");
            return result;
        }

        else if (mat_1.length == 1) { // for 1x1 matrix
            result = (int) mat_1[0][0];
            return result;
        }

        else if (mat_1.length == 2) { // for 2x2 matrix
            result = (((int) mat_1[0][0]) * ((int) mat_1[1][1]) - (((int) mat_1[0][1])) * ((int) mat_1[1][0]));
            return result;
        }

        else if (mat_1.length == 3) { // for 3x3 matrix
            for (int i = 0; i < mat_1[0].length; i++) {
                double[][] help = new double[mat_1.length - 1][mat_1[0].length - 1];

                for (int j = 1; j < mat_1.length; j++) {
                    for (int f = 0; f < mat_1[0].length; f++) {

                        if (f < i) {
                            help[j - 1][f] = mat_1[j][f];
                        } else if (f > i) {
                            help[j - 1][f - 1] = mat_1[j][f];
                        }
                    }
                }
                result += mat_1[0][i] * Math.pow(-1, (int) i) * determinant(help);
            }
        }
        return result;
    }

    public double[][] transpose(double[][] mat_1) {
        int rows = mat_1.length;
        int columns = mat_1[0].length;
        double[][] result = new double[columns][rows];
        for (int i = 0; i < rows; i++) {
            for (int j = 0; j < columns; j++) {
                result[j][i] = mat_1[i][j];
            }
        }
        return result;
    }

    public int determinant_2(double[][] mat_1) {
        int result = 0;
        int c = mat_1.length - 1;
        if (mat_1.length != mat_1[0].length) {
            System.out.println("error: matrix must be square matrix for determinant");
            return result;
        }

        else if (c == 1) { // for 1x1 matrix
            result = (int) mat_1[0][0];
            return result;
        }

        else if (c == 2) { // for 2x2 matrix
            result = (((int) mat_1[0][0]) * ((int) mat_1[1][1]) - (((int) mat_1[0][1])) * ((int) mat_1[1][0]));
            return result;
        }

        return result;
    }

    public double[][] cofactor(double[][] mat_1, double[][] help, int rows, int i, int j) {
        int a = 0, b = 0;
        for (int row = 0; row < rows; row++) {
            for (int col = 0; col < rows; col++) {
                if (row != i && col != j) {
                    help[a][b++] = mat_1[row][col];

                    if (b == rows - 1) {
                        b = 0;
                        a++;
                    }
                }
            }
        }
        return help;
    }

    public double[][] adjoint(double[][] mat_1) {
        int rows = mat_1.length;
        double[][] result = new double[rows][rows];

        if (rows == 1) { // for 1x1 matrix
            result[0][0] = 1;
        }

        double[][] help = new double[rows][rows];
        int sign = 1;

        for (int i = 0; i < rows; i++) {
            for (int j = 0; j < rows; j++) {

                help = cofactor(mat_1, help, rows, i, j);

                result[j][i] = (sign) * (determinant_2(help));
            }
        }

        for (int i = 0; i < rows; i++) {
            for (int j = 0; j < rows; j++) {
                if ((i + j) % 2 != 0) {
                    result[j][i] = (-1) * (result[j][i]);
                }
            }
        }

        return result;
    }

    public double[][] inverse(double[][] mat_1) {
        int rows = mat_1.length;
        double[][] result = new double[rows][rows];
        double[][] mat_1_adj = new double[rows][rows];

        int det = determinant(mat_1);
        if (det == 0) {
            System.out.print("It is a singular matrix and it has no inverse");
            result[0][0] = 0;
            return result;
        }
        mat_1_adj = adjoint(mat_1);
        for (int i = 0; i < rows; i++) {
            for (int j = 0; j < rows; j++) {
                result[i][j] = mat_1_adj[i][j] / (float) det;
            }
        }
        return result;
    }

    public void row_wise_mean(double[][] mat_1) {
        int rows = mat_1.length;
        int columns = mat_1[0].length;
        for (int i = 0; i < rows; i++) {
            int sum = 0;
            for (int j = 0; j < columns; j++) {
                sum += mat_1[i][j];
            }
            System.out.println("Row " + (i + 1) + " mean = " + (double) sum / columns);
        }
    }

    public void column_wise_mean(double[][] mat_1) {
        int rows = mat_1.length;
        int columns = mat_1[0].length;
        for (int i = 0; i < columns; i++) {
            int sum = 0;
            for (int j = 0; j < rows; j++) {
                sum += mat_1[j][i];
            }
            System.out.println("Column " + (i + 1) + " mean = " + (double) sum / rows);
        }
    }

    public void mean_of_all_elements(double[][] mat_1) {
        int rows = mat_1.length;
        int columns = mat_1[0].length;
        int sum = 0;
        for (int i = 0; i < rows; i++) {
            for (int j = 0; j < columns; j++) {
                sum += mat_1[i][j];
            }
        }
        System.out.println("mean of all elements = " + (double) sum / (rows * columns));
    }

    public double[][] singleton_addition(double[][] mat_1, double[][] mat_2) {
        int rows = mat_1.length;
        int columns = mat_1[0].length;
        float sing = (float) mat_2[0][0];
        double[][] result = new double[rows][columns];
        for (int i = 0; i < rows; i++) {
            for (int j = 0; j < columns; j++) {
                result[i][j] = mat_1[i][j] + sing;
            }
        }
        matrix_print(result);
        return result;
    }

    public double[][] singleton_subtraction(double[][] mat_1, double[][] mat_2) {
        int rows = mat_1.length;
        int columns = mat_1[0].length;
        float sing = (float) mat_2[0][0];
        double[][] result = new double[rows][columns];
        for (int i = 0; i < rows; i++) {
            for (int j = 0; j < columns; j++) {
                result[i][j] = mat_1[i][j] - sing;
            }
        }
        matrix_print(result);
        return result;
    }

    public double[][] singleton_multiplication(double[][] mat_1, double[][] mat_2) {
        int rows = mat_1.length;
        int columns = mat_1[0].length;
        double sing = mat_2[0][0];
        double[][] result = new double[rows][columns];
        for (int i = 0; i < rows; i++) {
            for (int j = 0; j < columns; j++) {
                result[i][j] = mat_1[i][j] * sing;
            }
        }
        matrix_print(result);
        return result;
    }

    public void singleton_division(double[][] mat_1, double[][] mat_2) {
        int rows = mat_1.length;
        int columns = mat_1[0].length;
        if ((int) mat_2[0][0] == 0) {
            System.out.println("error: a matrix can not be divisible by a singular matrix.");
            return;
        }
        float sing = (float) mat_2[0][0];
        double[][] result = new double[rows][columns];
        for (int i = 0; i < rows; i++) {
            for (int j = 0; j < columns; j++) {
                result[i][j] = mat_1[i][j] / sing;
            }
        }
        matrix_print(result);
    }

    public double[][] Compute_A_plus_A_transpose(double[][] mat_1) {
        int rows = mat_1.length;
        int columns = mat_1[0].length;
        double[][] result = new double[rows][columns];
        double[][] mat_1_transpose = new double[columns][rows];
        mat_1_transpose = transpose(mat_1);
        result = addition(mat_1, mat_1_transpose);
        return result;
    }

    public void solve_linear_equation(double[][] mat_1, double[][] mat_2) {
        if (mat_1.length != mat_1[0].length) {
            System.out.println("error: this matrix should be square matrix.");
            return;
        }

        if (mat_1.length != mat_2.length) {
            System.out.println("error: this matrix should have same no. of rows as previous one.");
            return;
        }

        double[][] result = new double[mat_1.length][1];
        double[][] mat_1_inverse = new double[mat_1.length][mat_1[0].length];
        mat_1_inverse = inverse(mat_1);
        result = multiplication(mat_1_inverse, mat_2);
    }

}

class matrix_type implements common_methods {

    public int determinant(double[][] mat_1) {
        int result = 0;

        if (mat_1.length == 1) { // for 1x1 matrix
            result = (int) mat_1[0][0];
            return result;
        }

        else if (mat_1.length == 2) { // for 2x2 matrix
            result = (((int) mat_1[0][0]) * ((int) mat_1[1][1]) - (((int) mat_1[0][1])) * ((int) mat_1[1][0]));
            return result;
        }

        else if (mat_1.length == 3) { // for 3x3 matrix
            for (int i = 0; i < mat_1[0].length; i++) {
                double[][] help = new double[mat_1.length - 1][mat_1[0].length - 1];

                for (int j = 1; j < mat_1.length; j++) {
                    for (int f = 0; f < mat_1[0].length; f++) {

                        if (f < i) {
                            help[j - 1][f] = mat_1[j][f];
                        } else if (f > i) {
                            help[j - 1][f - 1] = mat_1[j][f];
                        }
                    }
                }
                result += mat_1[0][i] * Math.pow(-1, (int) i) * determinant(help);
            }
        }
        return result;
    }

    public double[][] transpose(double[][] mat_1) {
        int rows = mat_1.length;
        int columns = mat_1[0].length;
        double[][] result = new double[columns][rows];
        for (int i = 0; i < rows; i++) {
            for (int j = 0; j < columns; j++) {
                result[j][i] = mat_1[i][j];
            }
        }
        return result;
    }

    public ArrayList<String> type_find(double[][] mat_type) {
        ArrayList<String> store_type_help = new ArrayList<String>();

        int rows = mat_type.length;
        int columns = mat_type[0].length;

        if (rows != columns) {
            store_type_help.add("Rectangular matrix");
        }

        if (rows == 1 && columns > 1) {
            store_type_help.add("Row matrix");
        }

        if (columns == 1 && rows > 1) {
            store_type_help.add("Column matrix");
        }

        if ((rows == columns) || (rows != 1 && columns != 1)) {
            store_type_help.add("Square matrix");
        }

        if (rows == columns) {
            int s = 0;
            for (int i = 0; i < rows; i++) {
                for (int j = 0; j < columns; j++) {
                    if ((int) mat_type[i][j] != (int) mat_type[j][i]) {
                        s = 1;
                        break;
                    }
                }
            }
            if (s == 0) {
                store_type_help.add("Symmetric matrix");
            }
        }

        if (rows == columns) {
            int s = 0;
            for (int i = 0; i < rows; i++) {
                for (int j = 0; j < columns; j++) {
                    if ((int) mat_type[i][j] != (-1) * (int) mat_type[j][i]) {
                        s = 1;
                        break;
                    }
                }
            }
            if (s == 0) {
                store_type_help.add("Skew symmetric matrix");
            }
        }

        if (rows == columns) {
            int s = 0;
            for (int i = 0; i < rows; i++) {
                for (int j = 0; j < columns; j++) {
                    if (i > j) {
                        if ((int) mat_type[i][j] != 0) {
                            s = 1;
                        }
                    }
                }
            }
            if (s == 0) {
                store_type_help.add("Upper triangular matrix");
            }
        }

        if (rows == columns) {
            int s = 0;
            for (int i = 0; i < rows; i++) {
                for (int j = 0; j < columns; j++) {
                    if (i < j) {
                        if ((int) mat_type[i][j] != 0) {
                            s = 1;
                        }
                    }
                }
            }
            if (s == 0) {
                store_type_help.add("Lower triangular matrix");
            }
        }

        if ((rows == columns) && ((int) determinant(mat_type) == 0)) {
            store_type_help.add("Singular matrix");
        }

        if (rows == columns) {
            int s = 0;
            for (int i = 0; i < rows; i++) {
                for (int j = 0; j < columns; j++) {
                    if (i != j) {
                        if ((int) mat_type[i][j] != 0) {
                            s = 1;
                        }
                    }
                }
            }
            if (s == 0) {
                store_type_help.add("Diagonal matrix");
            }
        }

        if (rows == columns) {
            int s = 0;
            int check = (int) mat_type[0][0];
            for (int i = 0; i < rows; i++) {
                for (int j = 0; j < columns; j++) {
                    if (i != j) {
                        if ((int) mat_type[i][j] != 0) {
                            s = 1;
                        }
                    }
                    if (i == j) {
                        if ((int) mat_type[i][j] != check) {
                            s = 1;
                        }
                    }
                }
            }
            if (s == 0) {
                store_type_help.add("Scalar matrix");
            }
        }

        if ((rows == columns) && ((int) determinant(mat_type) == 1)) {
            store_type_help.add("Identity matrix");
        }

        if (rows == 1 && columns == 1) {
            store_type_help.add("Singleton matrix");
        }

        if (rows == 1 && columns == 1) {
            int s = 0;
            for (int i = 0; i < rows; i++) {
                for (int j = 0; j < columns; j++) {
                    if ((int) mat_type[i][j] != 1) {
                        s = 1;
                    }
                }
            }
            if (s == 0) {
                store_type_help.add("Ones matrix");
            }
        }

        if (true) {
            int s = 0;
            for (int i = 0; i < rows; i++) {
                for (int j = 0; j < columns; j++) {
                    if ((int) mat_type[i][j] != 0) {
                        s = 1;
                    }
                }
            }
            if (s == 0) {
                store_type_help.add("Null matrix");
            }
        }

        return store_type_help;
    }

    static double get_random() {
        Random ran = new Random();
        int num = ran.nextInt(9);
        while (num == 0) {
            num = ran.nextInt(9);
        }
        double d = num;
        return d;
    }

    static double[][] default_matrix(int rows, int columns) {
        matrix mt = new matrix();

        double[][] result = new double[rows][columns];
        for (int i = 0; i < rows; i++) {
            for (int j = 0; j < columns; j++) {
                result[i][j] = get_random();
            }
        }
        mt.matrix_print(result);
        System.out.println();
        return result;
    }

    public double[][] default_matrix_generator(int given_type) {
        matrix mt = new matrix();

        if (given_type == 1) {
            double[][] result = new double[2][3];
            result = default_matrix(2, 3);
            return result;
        }

        if (given_type == 2) {
            double[][] result = new double[1][3];
            result = default_matrix(1, 3);
            return result;
        }

        if (given_type == 3) {
            double[][] result = new double[3][1];
            result = default_matrix(3, 1);
            return result;
        }

        if (given_type == 4) {
            double[][] result = new double[3][3];
            result = default_matrix(3, 3);
            return result;
        }

        if (given_type == 5) {
            int rows = 3;
            int columns = 3;
            double[][] result = new double[rows][columns];
            for (int i = 0; i < rows; i++) {
                for (int j = 0; j < columns; j++) {
                    if (i == j) {
                        result[i][j] = get_random();
                    }
                }
            }
            result[1][0] = get_random();
            result[2][0] = get_random();
            result[2][1] = get_random();

            result[0][1] = result[1][0];
            result[0][2] = result[2][0];
            result[1][2] = result[2][1];

            mt.matrix_print(result);
            System.out.println();
            return result;
        }

        if (given_type == 6) {
            int rows = 3;
            int columns = 3;
            double[][] result = new double[rows][columns];
            for (int i = 0; i < rows; i++) {
                for (int j = 0; j < columns; j++) {
                    if (i == j) {
                        result[i][j] = get_random();
                    }
                }
            }
            result[1][0] = get_random();
            result[2][0] = get_random();
            result[2][1] = get_random();

            result[0][1] = (-1) * result[1][0];
            result[0][2] = (-1) * result[2][0];
            result[1][2] = (-1) * result[2][1];

            mt.matrix_print(result);
            System.out.println();
            return result;
        }

        if (given_type == 7) {
            int rows = 3;
            int columns = 3;
            double[][] result = new double[rows][columns];
            for (int i = 0; i < rows; i++) {
                for (int j = 0; j < columns; j++) {
                    if (i > j) {
                        result[i][j] = 0;
                    } else {
                        result[i][j] = get_random();
                    }
                }
            }
            mt.matrix_print(result);
            System.out.println();
            return result;
        }

        if (given_type == 8) {
            int rows = 3;
            int columns = 3;
            double[][] result = new double[rows][columns];
            for (int i = 0; i < rows; i++) {
                for (int j = 0; j < columns; j++) {
                    if (i < j) {
                        result[i][j] = 0;
                    } else {
                        result[i][j] = get_random();
                    }
                }
            }
            mt.matrix_print(result);
            System.out.println();
            return result;
        }

        if (given_type == 9) {
            int rows = 2;
            int columns = 2;
            double[][] result = new double[rows][columns];
            result[0][0] = 3.0;
            result[0][1] = 6.0;
            result[1][0] = 2.0;
            result[1][1] = 4.0;

            mt.matrix_print(result);
            System.out.println();
            return result;
        }

        if (given_type == 10) {
            int rows = 3;
            int columns = 3;
            double[][] result = new double[rows][columns];
            for (int i = 0; i < rows; i++) {
                for (int j = 0; j < columns; j++) {
                    if (i == j) {
                        result[i][j] = get_random();
                    } else {
                        result[i][j] = 0;
                    }
                }
            }
            mt.matrix_print(result);
            System.out.println();
            return result;
        }

        if (given_type == 11) {
            int rows = 3;
            int columns = 3;
            double ele = get_random();
            double[][] result = new double[rows][columns];
            for (int i = 0; i < rows; i++) {
                for (int j = 0; j < columns; j++) {
                    if (i == j) {
                        result[i][j] = ele;
                    } else {
                        result[i][j] = 0;
                    }
                }
            }
            mt.matrix_print(result);
            System.out.println();
            return result;
        }

        if (given_type == 12) {
            int rows = 3;
            int columns = 3;
            double[][] result = new double[rows][columns];
            for (int i = 0; i < rows; i++) {
                for (int j = 0; j < columns; j++) {
                    if (i == j) {
                        result[i][j] = 1.0;
                    } else {
                        result[i][j] = 0;
                    }
                }
            }
            mt.matrix_print(result);
            System.out.println();
            return result;
        }

        if (given_type == 13) {
            int rows = 1;
            int columns = 1;
            double[][] result = new double[rows][columns];
            result[0][0] = get_random();
            mt.matrix_print(result);
            System.out.println();
            return result;
        }

        if (given_type == 14) {
            int rows = 3;
            int columns = 3;
            double[][] result = new double[rows][columns];
            for (int i = 0; i < rows; i++) {
                for (int j = 0; j < columns; j++) {
                    result[i][j] = 1.0;
                }
            }
            mt.matrix_print(result);
            System.out.println();
            return result;
        }

        if (given_type == 15) {
            int rows = 3;
            int columns = 3;
            double[][] result = new double[rows][columns];
            for (int i = 0; i < rows; i++) {
                for (int j = 0; j < columns; j++) {
                    result[i][j] = 0;
                }
            }
            mt.matrix_print(result);
            System.out.println();
            return result;
        } else {
            double[][] result = new double[1][1];
            return result;
        }
    }

}

public class MatrixCalculator {
    static Scanner scn = new Scanner(System.in);

    static void task_list() {
        System.out.println("\n-------------------------------------------------------------");
        System.out.println("\n               TASK LIST              ");
        System.out.println("-------------------------------------------------------------");
        System.out.println("1. Create matrix.\r\n"
                + "2. Create matrices of requested matrix-types.\r\n"
                + "3. Change the elements of a matrix.\r\n"
                + "4. Display all the matrix-type of a requested matrix.\r\n"
                + "5. Perform addition, subtraction, multiplication & division.\r\n"
                + "6. Perform element-wise operations.\r\n"
                + "7. Transpose matrices.\r\n"
                + "8. Inverse matrices.\r\n"
                + "9. Compute means.\r\n\""
                + "10. Compute determinants.\r\n"
                + "11. Use singleton matrices as scalars, if requested.\r\n"
                + "12. Compute A+A^T for a matrix A.\r\n"
                + "13. Compute Eigen vectors and values.\r\n"
                + "14. Solve sets of linear equations using matrices.\r\n"
                + "15. Retrieve all the existing matrices (entered or created) having requested matrix-type labels.\r\n"
                + "16. Exit");
        System.out.println("-------------------------------------------------------------");
    }

    public static void main(String[] args) {

        matrix mat = new matrix();
        matrix_type tye = new matrix_type();
        ArrayList<double[][]> store_matrix = new ArrayList<double[][]>();
        ArrayList<ArrayList<String>> store_type = new ArrayList<ArrayList<String>>();

        while (true) {
            System.out.println(" ");
            System.out.println("select the task you want to perform: ");
            task_list();

            int task = scn.nextInt();

            if (task == 1) {
                System.out.print("enter no. of rows: ");
                int row = scn.nextInt();
                System.out.print("enter no. of columns: ");
                int col = scn.nextInt();
                System.out.println();

                store_matrix.add(mat.matrix_generator(row, col));

                store_type.add(tye.type_find(store_matrix.get(store_matrix.size() - 1)));
            }

            else if (task == 2) {
                System.out.println("select the matrix type: ");
                System.out.print("1. Rectangular Matrix\r\n"
                        + "2. Row Matrix\r\n"
                        + "3. Column Matrix\r\n"
                        + "4. Square Matrix\r\n"
                        + "5. Symmetric Matrix\r\n"
                        + "6. Skew-symmetric Matrix\r\n"
                        + "7. Upper-triangular Matrix\r\n"
                        + "8. Lower-triangular Matrix\r\n"
                        + "9. Singular Matrix\r\n"
                        + "10. Diagonal Matrix\r\n"
                        + "11. Scalar Matrix\r\n"
                        + "12. Identity Matrix\r\n"
                        + "13. Singleton Matrix\r\n"
                        + "14. Ones Matrix\r\n"
                        + "15. Null Matrix");

                System.out.println(" ");
                int select = scn.nextInt();
                store_matrix.add(tye.default_matrix_generator(select));

                store_type.add(tye.type_find(store_matrix.get(store_matrix.size() - 1)));
            }

            else if (task == 3) {
                System.out.println("select a matrix (number-wise) to change its required elements: ");
                int m1 = scn.nextInt();
                System.out.println("Tell the no. of elements you want to change: ");
                int number = scn.nextInt();

                int rows = store_matrix.get(m1 - 1).length;
                int columns = store_matrix.get(m1 - 1)[0].length;

                double[][] copy = new double[rows][columns];
                copy = store_matrix.get(m1 - 1);
                ArrayList<String> check = new ArrayList<String>();

                for (int k = 0; k < number; k++) {
                    System.out.println("select one of the dimension of the selected matrix to change: ");
                    System.out.println("select its row dimension: ");
                    int r = scn.nextInt();
                    System.out.println("select its column dimension: ");
                    int c = scn.nextInt();

                    for (int i = 0; i < rows; i++) {
                        for (int j = 0; j < columns; j++) {
                            if ((i == (r - 1)) && (j == (c - 1))) {
                                System.out.println("tell new element value: ");
                                int ask = scn.nextInt();
                                double ask_2 = ask;
                                copy[i][j] = ask_2;
                            }
                        }
                    }
                }
                check = tye.type_find(copy);

                int s = 0;
                for (int i = 0; i < store_type.get(m1 - 1).size(); i++) {
                    for (int j = 0; j < check.size(); j++) {
                        if (store_type.get(m1 - 1).get(i) == check.get(j)) {
                            s++;
                        }
                    }
                }

                if (s == store_type.get(m1 - 1).size()) {
                    store_type.set((m1 - 1), check);
                    store_matrix.set((m1 - 1), copy);
                    mat.matrix_print(copy);
                } else {
                    System.out.println(
                            "we can not change this element by the given value because it changes the matrix type. ");
                    System.out.println("");
                }

            }

            else if (task == 4) {
                System.out.println("which matrix do you want to select to see its types. ");
                System.out.println("select number-wise matrix: ");
                int m1 = scn.nextInt();
                System.out.println("Type of required matrix: ");
                for (int i = 0; i < store_type.get(m1 - 1).size(); i++) {
                    System.out.println(store_type.get(m1 - 1).get(i));
                }
            }

            else if (task == 5) {
                System.out.println("which one do you want to perform: ");
                System.out.println("1. addition\r\n"
                        + "2. subtraction\r\n"
                        + "3. multiplication\r\n"
                        + "4. division");
                int select = scn.nextInt();

                System.out.println("select 2 matrices (number-wise): ");
                int m1 = scn.nextInt();
                int m2 = scn.nextInt();

                if (select == 1) {
                    mat.addition(store_matrix.get(m1 - 1), store_matrix.get(m2 - 1));
                } else if (select == 2) {
                    mat.subtraction(store_matrix.get(m1 - 1), store_matrix.get(m2 - 1));
                } else if (select == 3) {
                    mat.multiplication(store_matrix.get(m1 - 1), store_matrix.get(m2 - 1));
                } else if (select == 4) {
                    mat.division(store_matrix.get(m1 - 1), store_matrix.get(m2 - 1));
                }

            }

            else if (task == 6) {
                System.out.println("which element wise operation do you want to perform: ");
                System.out.println("1. addition\r\n"
                        + "2. subtraction\r\n"
                        + "3. multiplication\r\n"
                        + "4. division");
                int select = scn.nextInt();

                System.out.println("select 2 matrices (number-wise): ");
                int m1 = scn.nextInt();
                int m2 = scn.nextInt();

                if (select == 1) {
                    mat.addition(store_matrix.get(m1 - 1), store_matrix.get(m2 - 1));
                } else if (select == 2) {
                    mat.subtraction(store_matrix.get(m1 - 1), store_matrix.get(m2 - 1));
                } else if (select == 3) {
                    mat.element_wise_multiplication(store_matrix.get(m1 - 1), store_matrix.get(m2 - 1));
                } else if (select == 4) {
                    mat.element_wise_division(store_matrix.get(m1 - 1), store_matrix.get(m2 - 1));
                }
            }

            else if (task == 7) {
                System.out.println("select a matrix (number-wise) to transpose: ");
                int m1 = scn.nextInt();

                double[][] help = new double[store_matrix.get(m1 - 1).length][store_matrix.get(m1 - 1)[0].length];

                help = mat.transpose(store_matrix.get(m1 - 1));
                mat.matrix_print(help);
            }

            else if (task == 8) {
                System.out.println("select a matrix (number-wise) to inverse: ");
                int m1 = scn.nextInt();

                double[][] help = new double[store_matrix.get(m1 - 1).length][store_matrix.get(m1 - 1)[0].length];

                help = mat.inverse(store_matrix.get(m1 - 1));
                mat.matrix_print(help);
            }

            else if (task == 9) {
                System.out.println("which type of mean do you want to perform: ");
                System.out.println("1. row-wise mean\r\n"
                        + "2. column-wise mean\r\n"
                        + "3. mean of all the elements");
                int select = scn.nextInt();

                System.out.println("select a matrix (number-wise) to perform mean: ");
                int m1 = scn.nextInt();

                if (select == 1) {
                    mat.row_wise_mean(store_matrix.get(m1 - 1));
                } else if (select == 2) {
                    mat.column_wise_mean(store_matrix.get(m1 - 1));
                } else if (select == 3) {
                    mat.mean_of_all_elements(store_matrix.get(m1 - 1));
                }
            }

            else if (task == 10) {
                System.out.println("select a matrix (number-wise) to calculate determinant: ");
                int m1 = scn.nextInt();

                int help = mat.determinant(store_matrix.get(m1 - 1));
                System.out.println("determinant is: " + help);
            }

            else if (task == 11) {
                System.out.println("Do you allow using singleton matrices as a scalar value?");
                System.out.println("1. yes\r\n"
                        + "2. no");
                int option = scn.nextInt();

                if (option == 1) {
                    System.out.println("which operation do you want to perform using singleton matrix: ");
                    System.out.println("1. addition\r\n"
                            + "2. subtraction\r\n"
                            + "3. multiplication\r\n"
                            + "4. division");
                    int select = scn.nextInt();

                    System.out.println("select 2 matrices (number-wise): ");
                    System.out.println("non-singleton matrix: ");
                    int m1 = scn.nextInt();
                    System.out.println("singleton matrix: ");
                    int m2 = scn.nextInt();

                    if (select == 1) {
                        mat.singleton_addition(store_matrix.get(m1 - 1), store_matrix.get(m2 - 1));
                    } else if (select == 2) {
                        mat.singleton_subtraction(store_matrix.get(m1 - 1), store_matrix.get(m2 - 1));
                    } else if (select == 3) {
                        mat.singleton_multiplication(store_matrix.get(m1 - 1), store_matrix.get(m2 - 1));
                    } else if (select == 4) {
                        mat.singleton_division(store_matrix.get(m1 - 1), store_matrix.get(m2 - 1));
                    }
                }
            }

            else if (task == 12) {
                System.out.println("select a matrix (number-wise) to calculate A+A': ");
                int m1 = scn.nextInt();

                mat.Compute_A_plus_A_transpose(store_matrix.get(m1 - 1));
            }

            else if (task == 13) {
                System.out.println("select a matrix (number-wise) to find eigen values: ");
                int m1 = scn.nextInt();

                double[][] help = new double[store_matrix.get(m1 - 1).length][store_matrix.get(m1 - 1)[0].length];

                mat.eigen_values(store_matrix.get(m1 - 1));
            }

            else if (task == 14) {
                System.out.println("Solving a set of linear equation. ");
                System.out.println("select a square matrix (number-wise): ");
                int m1 = scn.nextInt();
                System.out.println("select a column matrix (number-wise): ");
                int m2 = scn.nextInt();

                mat.solve_linear_equation(store_matrix.get(m1 - 1), store_matrix.get(m2 - 1));
            }

            else if (task == 15) {
                System.out.println("select a type of matrix to rerieve all matrices.");
                System.out.println("name should be in format, eg. Row matrix.");
                String name = scn.nextLine();
                if ((name.equals(null)) || (name.equals(""))) {
                    name = scn.nextLine();
                }

                for (int i = 0; i < store_type.size(); i++) {
                    for (int j = 0; j < store_type.get(i).size(); j++) {
                        if (name == store_type.get(i).get(j)) {
                            System.out.println("matrix: ");
                            mat.matrix_print(store_matrix.get(i));
                            System.out.println(" ");
                        }
                    }
                }

            }

            else if (task == 16) {
                break;
            }

        }

    }
}