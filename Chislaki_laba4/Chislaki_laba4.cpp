#include <iostream>
#include <iomanip>
using namespace std;

double** input_Table(double** Table_Matrix, int amount_measures)
{
    cout << "Введите таблицу измерений: " << endl;
	for (int i = 0; i < 2; i++)
	{
		Table_Matrix[i] = new double [amount_measures];
		for (int j = 0; j < amount_measures; j++)
		{
			cout << "Table_Matrix[" << i << "][" << j << "]= ";
			cin >> Table_Matrix[i][j];
		}
	}
	return Table_Matrix;
	for (int i = 0; i < 2; i++)
	{
		delete[]Table_Matrix[i];
	}
}
double* Solution_by_the_Gaussian_method(double** Matrix_A, double* Vector_b, int rank)
{
    double* Vector_X, max;
    int k, index;
    const double eps = 0.000001;  // точность
    Vector_X = new double[rank];
    k = 0;
    while (k < rank)
    {
        // Поиск строки с максимальным a[i][k]
        max = abs(Matrix_A[k][k]);
        index = k;
        for (int i = k + 1; i < rank; i++)
        {
            if (abs(Matrix_A[i][k]) > max)
            {
                max = abs(Matrix_A[i][k]);
                index = i;
            }
        }
        if (max < eps)
        {
            // нет ненулевых диагональных элементов
            cout << "Решение получить невозможно из-за нулевого столбца ";
            cout << index << " матрицы A" << endl;
            return 0;
        }
        // Перестановка строк
        for (int j = 0; j < rank; j++)
        {
            double temp = Matrix_A[k][j];
            Matrix_A[k][j] = Matrix_A[index][j];
            Matrix_A[index][j] = temp;
        }
        double temp = Vector_b[k];
        Vector_b[k] = Vector_b[index];
        Vector_b[index] = temp;
        // Нормализация уравнений
        for (int i = k; i < rank; i++)
        {
            double temp = Matrix_A[i][k];
            if (abs(temp) < eps) continue; // для нулевого коэффициента пропустить
            for (int j = 0; j < rank; j++)
                Matrix_A[i][j] = Matrix_A[i][j] / temp;
            Vector_b[i] = Vector_b[i] / temp;
            if (i == k)  continue; // уравнение не вычитать само из себя
            for (int j = 0; j < rank; j++)
                Matrix_A[i][j] = Matrix_A[i][j] - Matrix_A[k][j];
            Vector_b[i] = Vector_b[i] - Vector_b[k];
        }
        k++;
    }
    // обратная подстановка
    for (k = rank - 1; k >= 0; k--)
    {
        Vector_X[k] = Vector_b[k];
        for (int i = 0; i < k; i++)
            Vector_b[i] = Vector_b[i] - Matrix_A[i][k] * Vector_X[k];
    }
    return Vector_X;
    delete[]Vector_X;
}
void output_Table(double** Table_Matrix, double* Vector_X, int amount_measures)
{
    double* Result;
    Result = new double[amount_measures];
    for (int i = 0; i < amount_measures; i++)
    {
        Result[i] = Vector_X[0] / Table_Matrix[0][i] + Vector_X[1];
    }
    cout << "Напор жидкости H:    Коэффициент истечения mu(табличный):   Коэффициент истечения mu(вычисленный):" << endl;
    for (int i = 0; i < amount_measures; i++)
    {
        cout << setw(17) << Table_Matrix[0][i] << setw(40) << Table_Matrix[1][i] <<setw(41) << Result[i] << endl;
    }
    delete[]Result;
}
int main()
{
	setlocale(LC_ALL, "rus");

	int amount_measures;//6
	cout << "Введите количество измерений: " << endl;
	cin >> amount_measures;

	double** Table_matrix;
	Table_matrix = new double* [amount_measures];
	Table_matrix = input_Table(Table_matrix, amount_measures);
    // замена переменных для линеаризации
    //////////////////////////////////////////////////////////////
    double* Array_Xi;// массив иксов
    Array_Xi = new double[amount_measures];
    for (int i = 0; i < amount_measures; i++)
    {
        Array_Xi[i] = 1 / Table_matrix[0][i];
    }
    /////////////////////////////////////////////////////////////
    double* Array_Xi_squares;// массив квдратов иксов
    Array_Xi_squares = new double[amount_measures];
    for (int i = 0; i < amount_measures; i++)
    {
        Array_Xi_squares[i] = Array_Xi[i] * Array_Xi[i];
    }
    //////////////////////////////////////////////////////////////////
    double* Array_Yi;// массив игриков
    Array_Yi = new double[amount_measures];
    for (int i = 0; i < amount_measures; i++)
    {
        Array_Yi[i] = Table_matrix[1][i];
    }
    ///////////////////////////////////////////////////////////////////
    double* Array_XiYi;// массив произведений иксов и игриков
    Array_XiYi = new double[amount_measures];
    for (int i = 0; i < amount_measures; i++)
    {
        Array_XiYi[i] = Array_Xi[i] * Array_Yi[i];
    }
    //////////////////////////////////////////////////////////////////////

    double Sum_squares_Xi = 0;
    for (int i = 0; i < amount_measures; i++)
    {
        Sum_squares_Xi += Array_Xi_squares[i];
    }

    double Sum_Xi = 0;
    for (int i = 0; i < amount_measures; i++)
    {
        Sum_Xi += Array_Xi[i];
    }

    double Sum_XiYi = 0;
    for (int i = 0; i < amount_measures; i++)
    {
        Sum_XiYi += Array_XiYi[i];
    }

    double Sum_Yi = 0;
    for (int i = 0; i < amount_measures; i++)
    {
        Sum_Yi += Array_Yi[i];
    }
    ////////////////////////////////////////////////////////////

    double A, B;
    double** Matrix_A, * Vector_b, * Vector_X;// для решения СЛАУ методом Гаусса

    Matrix_A = new double* [2];
    for (int i = 0; i < 1; i++)
    {
        Matrix_A[i] = new double[2];
        Matrix_A[i][0] = Sum_squares_Xi;
        Matrix_A[i][1] = Sum_Xi;
    }
    for (int i = 1; i < 2; i++)
    {
        Matrix_A[i] = new double[2];
        Matrix_A[i][0] = Sum_Xi;
        Matrix_A[i][1] = amount_measures;
    }
   /* Matrix_A[0] = new double[2];
    Matrix_A[0][0] = Sum_squares_Xi;
    Matrix_A[0][1] = Sum_Xi;
    Matrix_A[1] = new double[2];
    Matrix_A[1][0] = Sum_Xi;
    Matrix_A[1][1] = amount_measures;*/

    Vector_b = new double[2];
    Vector_b[0] = Sum_XiYi;
    Vector_b[1] = Sum_Yi;

    Vector_X = Solution_by_the_Gaussian_method(Matrix_A, Vector_b, 2);

    cout << "Ответ:" << endl
        << "A = " << Vector_X[0] << "  B = " << Vector_X[1] << endl;
    output_Table(Table_matrix, Vector_X, amount_measures);


    for (int i = 0; i < 2; i++)
    {
        delete[]Matrix_A[i];
    }
    delete[]Matrix_A;
    delete[]Vector_b;
    delete[]Array_Xi;
    delete[]Array_Yi;
    delete[]Array_Xi_squares;
    delete[]Array_XiYi;

	delete[]Table_matrix;
}