// carbon.cpp : このファイルには 'main' 関数が含まれています。プログラム実行の開始と終了がそこで行われます。
//

#include <stdio.h>
#include<math.h>

#define N 3 /* 逆行列を求める行列の行数・列数 */
#define MAX_ERR 1e-10 /* 許容する誤差 */
#ifndef M_PI
#define M_PI 3.14159265358979323846264338327950288
#endif


void Display(double Q[3][3][100], int I);
int check(double mat[N][N][100], double inv[N][N][100], int I);
void Syokika(double Q[3][3][100], int I);
double Tenchi(double Q[3][3]);
void GyakuGyoretsu(double Q[3][3][100], int I);
void Gyoretsu_Kakesan(double Q[3][3][100], double R[3][3][100], double result[3][3][100], int I);
void Input(double A[3][3][100], int I);
void Stress_Strain(double T_sigma[3][3][100], double Q[3][3][100], double T_epsilon[3][3][100], double result[3][3][100], int I);
void In_plane_stiffness_matrix_function(double Q_bar[3][3][100], double Thickness[100], double In_plane_stiffness_matrix[3][3][100], int I);//Iは(1からn-1まで)
void Tensile_bending_coupling_matrix_function(double Q_bar[3][3][100], double Thickness[100], double Tensile_bending_coupling_matrix[3][3][100], int I);//Iは(1からn-1まで)
void Out_of_plane_stiffness_matrix_function(double Q_bar[3][3][100], double Thickness[100], double Out_of_plane_stiffness_matrix[3][3][100], int I);//Iは(1からn-1まで)

int main()
{
	int n, i, s;
	double E_L, E_T, nu_LT, nu_TL, G_LT, Q[3][3][100], T_sigma[3][3][100], T_epsilon[3][3][100], theta[100], h[100], height_sum;
	double In_plane_stiffness_matrix[3][3][100], Tensile_bending_coupling_matrix[3][3][100], Out_of_plane_stiffness_matrix[3][3][100];

	Syokika(In_plane_stiffness_matrix, 99);
	Syokika(Tensile_bending_coupling_matrix, 99);
	Syokika(Out_of_plane_stiffness_matrix, 99);
	height_sum = 0;


	printf("何層の積層を扱うか？：");
	scanf_s("%d", &n);
	printf("疑似当方積層か？Yesなら2,Noなら1を代入せよ :");
	scanf_s("%d", &s);
	/*剛性行列[Q]を生成*/
	printf("E_L(繊維方向の弾性係数): ");//繊維方向の弾性係数
	scanf_s("%lf", &E_L);
	printf("E_T(繊維直行方向の弾性係数): ");//繊維直行方向の弾性係数
	scanf_s("%lf", &E_T);
	printf("nu_LT(主ポアソン比): ");//主ポアソン比
	scanf_s("%lf", &nu_LT);
	printf("nu_TL(従ポアソン比): ");//従ポアソン比
	scanf_s("%lf", &nu_TL);
	printf("G_LT(せん断弾性率): ");//せん断弾性率
	scanf_s("%lf", &G_LT);


	for (i = 0; i < n; i++)//層数分の剛性行列を生成する
	{
		int k;
		k = s * n - i - 1;
		if (s == 2)
		{
			Syokika(Q, i);
			Syokika(Q, k);
			Q[0][0][i] = Q[0][0][k] = E_L / (1 - (nu_LT * nu_TL));
			Q[0][1][i] = Q[0][1][k] = nu_LT * E_T / (1 - (nu_LT * nu_TL));
			Q[1][0][i] = Q[1][0][k] = nu_TL * E_L / (1 - (nu_LT * nu_TL));
			Q[1][1][i] = Q[1][1][k] = E_T / (1 - (nu_LT * nu_TL));
			Q[2][2][i] = Q[2][2][k] = G_LT;

			printf("%d層のh(積層板の厚さ)[mm]: ", i + 1);
			scanf_s("%lf", &h[i]);
		}
		else
		{
			Syokika(Q, i);
			Q[0][0][i] = E_L / (1 - (nu_LT * nu_TL));
			Q[0][1][i] = nu_LT * E_T / (1 - (nu_LT * nu_TL));
			Q[1][0][i] = nu_TL * E_L / (1 - (nu_LT * nu_TL));
			Q[1][1][i] = E_T / (1 - (nu_LT * nu_TL));
			Q[2][2][i] = G_LT;

			printf("%d層のh(積層板の厚さ)[mm]: ", i + 1);
			scanf_s("%lf", &h[i]);

		}
		if (s == 2)//対称積層の場合
		{
			h[k] = h[i];
			height_sum += h[i] * 2;
		}
		if (s == 1)//非対称積層の場合
		{
			height_sum += h[i];
		}
	}

	double height, h_r[100];
	if (s == 2)
	{
		for (i = 0; i < s * n + 1; i++)
		{
			if (i == 0)
			{
				h_r[i] = -height_sum / 2;
			}
			else
			{
				h_r[i] = h[i - 1] + h_r[i - 1];
			}
			//printf("h_r[%d]: %lf\n", i, h_r[i]);
		}

	}
	else
	{
		for (i = 0; i < n + 1; i++)
		{
			if (i == 0)
			{
				h_r[i] = -height_sum / 2;
			}
			else
			{
				h_r[i] = h[i - 1] + h_r[i - 1];
			}
		/*	printf("h_r[%d]: %lf\n", i, h_r[i]);*/
		}

	}

	for (i = 0; i < n; i++)
	{
		printf("Theta: ");
		scanf_s("%lf", &theta[i]);
		theta[i] = theta[i] * M_PI / 180;
		if (s == 2)
		{
			/*θだけ回転した座標系への応力の変換行列*/
			int k;
			k = s * n - i - 1;
			T_sigma[0][0][i] = T_sigma[0][0][k] = pow(cos(theta[i]), 2);
			T_sigma[0][1][i] = T_sigma[0][1][k] = pow(sin(theta[i]), 2);
			T_sigma[0][2][i] = T_sigma[0][2][k] = sin(2 * theta[i]);
			T_sigma[1][0][i] = T_sigma[1][0][k] = pow(sin(theta[i]), 2);
			T_sigma[1][1][i] = T_sigma[1][1][k] = pow(cos(theta[i]), 2);
			T_sigma[1][2][i] = T_sigma[1][2][k] = -sin(2 * theta[i]);
			T_sigma[2][0][i] = T_sigma[2][0][k] = -0.5 * sin(2 * theta[i]);
			T_sigma[2][1][i] = T_sigma[2][1][k] = 0.5 * sin(2 * theta[i]);
			T_sigma[2][2][i] = T_sigma[2][2][k] = cos(2 * theta[i]);

			/*θだけ回転した座標系へのひずみの変換行列*/
			T_epsilon[0][0][i] = T_epsilon[0][0][k] = pow(cos(theta[i]), 2);
			T_epsilon[0][1][i] = T_epsilon[0][1][k] = pow(sin(theta[i]), 2);
			T_epsilon[0][2][i] = T_epsilon[0][2][k] = 0.5 * sin(2 * theta[i]);
			T_epsilon[1][0][i] = T_epsilon[1][0][k] = pow(sin(theta[i]), 2);
			T_epsilon[1][1][i] = T_epsilon[1][1][k] = pow(cos(theta[i]), 2);
			T_epsilon[1][2][i] = T_epsilon[1][2][k] = -0.5 * sin(2 * theta[i]);
			T_epsilon[2][0][i] = T_epsilon[2][0][k] = -sin(2 * theta[i]);
			T_epsilon[2][1][i] = T_epsilon[2][1][k] = sin(2 * theta[i]);
			T_epsilon[2][2][i] = T_epsilon[2][2][k] = cos(2 * theta[i]);
		}
		else
		{
			/*θだけ回転した座標系への応力の変換行列*/
			T_sigma[0][0][i] = pow(cos(theta[i]), 2);
			T_sigma[0][1][i] = pow(sin(theta[i]), 2);
			T_sigma[0][2][i] = sin(2 * theta[i]);
			T_sigma[1][0][i] = pow(sin(theta[i]), 2);
			T_sigma[1][1][i] = pow(cos(theta[i]), 2);
			T_sigma[1][2][i] = -sin(2 * theta[i]);
			T_sigma[2][0][i] = -0.5 * sin(2 * theta[i]);
			T_sigma[2][1][i] = 0.5 * sin(2 * theta[i]);
			T_sigma[2][2][i] = cos(2 * theta[i]);

			/*θだけ回転した座標系へのひずみの変換行列*/
			T_epsilon[0][0][i] = pow(cos(theta[i]), 2);
			T_epsilon[0][1][i] = pow(sin(theta[i]), 2);
			T_epsilon[0][2][i] = 0.5 * sin(2 * theta[i]);
			T_epsilon[1][0][i] = pow(sin(theta[i]), 2);
			T_epsilon[1][1][i] = pow(cos(theta[i]), 2);
			T_epsilon[1][2][i] = -0.5 * sin(2 * theta[i]);
			T_epsilon[2][0][i] = -sin(2 * theta[i]);
			T_epsilon[2][1][i] = sin(2 * theta[i]);
			T_epsilon[2][2][i] = cos(2 * theta[i]);
		}

	}

	double Q_bar[3][3][100];
	if (s == 2)
	{
		for (i = 0; i < n * s; i++)
		{

			Syokika(In_plane_stiffness_matrix, i);
			Syokika(Tensile_bending_coupling_matrix, i);
			Syokika(Out_of_plane_stiffness_matrix, i);

			/*x-y座標系の応力とひずみの関係式の係数行列Q_barを求める*/

			Syokika(Q_bar, i);
			Stress_Strain(T_sigma, Q, T_epsilon, Q_bar, i);

			/*計算が正しいかチェック
			double Q_bar_check[3][3], l, m;
			Syokika(Q_bar_check);
			l = cos(theta);
			m = sin(theta);
			Q_bar_check[0][0] = Q[0][0] * pow(l, 4) + 2 * (Q[0][1] + 2 * Q[2][2]) * pow(l * m, 2) + Q[1][1] * pow(m, 4);
			Q_bar_check[1][0] = Q_bar_check[0][1] = Q[0][1] * (pow(l, 4) + pow(m, 4)) + (Q[0][0] + Q[1][1] - 4 * Q[2][2]) * pow(l * m, 2);
			Q_bar_check[1][1] = Q[0][0] * pow(m, 4) + 2 * (Q[0][1] + 2 * Q[2][2]) * pow(l * m, 2) + Q[1][1] * pow(l, 4);
			Q_bar_check[0][2] = Q_bar_check[2][0] = (Q[0][0] - Q[0][1] - 2 * Q[2][2]) * pow(l, 3) * m - (Q[1][1] - Q[0][1] - 2 * Q[2][2]) * l * pow(m, 3);
			Q_bar_check[1][2] = Q_bar_check[2][1] = (Q[0][0] - Q[0][1] - 2 * Q[2][2]) * l * pow(m, 3) - (Q[1][1] - Q[0][1] - 2 * Q[2][2]) * m * pow(l, 3);
			Q_bar_check[2][2] = (Q[0][0] + Q[1][1] - 2 * Q[0][1] - 2 * Q[2][2]) * pow(l * m, 2) + Q[2][2] * (pow(l, 4) + pow(m, 4));
			Display(Q_bar_check);
			*///ここまでの計算に間違いはない事を確認しているため、これ以後の計算に何らかの不具合が含まれている。


			/*printf("面内剛性行列[A]を表示\n=====================\n");
			Display(In_plane_stiffness_matrix);
			printf("=====================\n");
			printf("引っ張り-曲げカップリング行列[B]を表示\n=====================\n");
			Display(Tensile_bending_coupling_matrix);
			printf("=====================\n");
			printf("面外剛性行列[D]を表示\n=====================\n");
			Display(Out_of_plane_stiffness_matrix);
			printf("=====================\n");*/

		}
	}
	else
	{
		for (i = 0; i < n; i++)
		{

			Syokika(In_plane_stiffness_matrix, i);
			Syokika(Tensile_bending_coupling_matrix, i);
			Syokika(Out_of_plane_stiffness_matrix, i);

			/*x-y座標系の応力とひずみの関係式の係数行列Q_barを求める*/

			Syokika(Q_bar, i);
			Stress_Strain(T_sigma, Q, T_epsilon, Q_bar, i);

			/*計算が正しいかチェック
			double Q_bar_check[3][3], l, m;
			Syokika(Q_bar_check);
			l = cos(theta);
			m = sin(theta);
			Q_bar_check[0][0] = Q[0][0] * pow(l, 4) + 2 * (Q[0][1] + 2 * Q[2][2]) * pow(l * m, 2) + Q[1][1] * pow(m, 4);
			Q_bar_check[1][0] = Q_bar_check[0][1] = Q[0][1] * (pow(l, 4) + pow(m, 4)) + (Q[0][0] + Q[1][1] - 4 * Q[2][2]) * pow(l * m, 2);
			Q_bar_check[1][1] = Q[0][0] * pow(m, 4) + 2 * (Q[0][1] + 2 * Q[2][2]) * pow(l * m, 2) + Q[1][1] * pow(l, 4);
			Q_bar_check[0][2] = Q_bar_check[2][0] = (Q[0][0] - Q[0][1] - 2 * Q[2][2]) * pow(l, 3) * m - (Q[1][1] - Q[0][1] - 2 * Q[2][2]) * l * pow(m, 3);
			Q_bar_check[1][2] = Q_bar_check[2][1] = (Q[0][0] - Q[0][1] - 2 * Q[2][2]) * l * pow(m, 3) - (Q[1][1] - Q[0][1] - 2 * Q[2][2]) * m * pow(l, 3);
			Q_bar_check[2][2] = (Q[0][0] + Q[1][1] - 2 * Q[0][1] - 2 * Q[2][2]) * pow(l * m, 2) + Q[2][2] * (pow(l, 4) + pow(m, 4));
			Display(Q_bar_check);
			*///ここまでの計算に間違いはない事を確認しているため、これ以後の計算に何らかの不具合が含まれている。


			/*printf("面内剛性行列[A]を表示\n=====================\n");
			Display(In_plane_stiffness_matrix);
			printf("=====================\n");
			printf("引っ張り-曲げカップリング行列[B]を表示\n=====================\n");
			Display(Tensile_bending_coupling_matrix);
			printf("=====================\n");
			printf("面外剛性行列[D]を表示\n=====================\n");
			Display(Out_of_plane_stiffness_matrix);
			printf("=====================\n");*/

		}
	}


	if (s == 2)
	{
		for (i = 0; i < n * s; i++)
		{
			/*面内剛性行列(In_plane_stiffness_matrix)を求める*/
			In_plane_stiffness_matrix_function(Q_bar, h_r, In_plane_stiffness_matrix, i);

			/*引っ張り-曲げカップリング行列を求める*/
			Tensile_bending_coupling_matrix_function(Q_bar, h_r, Tensile_bending_coupling_matrix, i);

			/*面外剛性行列*/
			Out_of_plane_stiffness_matrix_function(Q_bar, h_r, Out_of_plane_stiffness_matrix, i);

		}
	}
	else
	{
		for (i = 0; i < n; i++)
		{
			/*面内剛性行列(In_plane_stiffness_matrix)を求める*/
			In_plane_stiffness_matrix_function(Q_bar, h_r, In_plane_stiffness_matrix, i);

			/*引っ張り-曲げカップリング行列を求める*/
			Tensile_bending_coupling_matrix_function(Q_bar, h_r, Tensile_bending_coupling_matrix, i);

			/*面外剛性行列*/
			Out_of_plane_stiffness_matrix_function(Q_bar, h_r, Out_of_plane_stiffness_matrix, i);

		}
	}


	printf("面内剛性行列[A]を表示\n=====================\n");
	Display(In_plane_stiffness_matrix, 99);
	printf("=====================\n");
	printf("引っ張り-曲げカップリング行列[B]を表示\n=====================\n");
	Display(Tensile_bending_coupling_matrix, 99);
	printf("=====================\n");
	printf("面外剛性行列[D]を表示\n=====================\n");
	Display(Out_of_plane_stiffness_matrix, 99);
	printf("=====================\n");


	GyakuGyoretsu(In_plane_stiffness_matrix, 99);
	printf("面内弾性特性: \n");
	printf("x方向のヤング率: %lf GPa\n", pow(In_plane_stiffness_matrix[0][0][99], -1) / height_sum);
	printf("y方向のヤング率: %lf GPa\n", pow(In_plane_stiffness_matrix[1][1][99], -1) / height_sum);
	printf("x方向のポアソン比: %lf\n", -In_plane_stiffness_matrix[1][0][99] * pow(In_plane_stiffness_matrix[0][0][99], -1));
	printf("y方向のポアソン比: %lf\n", -In_plane_stiffness_matrix[0][1][99] * pow(In_plane_stiffness_matrix[1][1][99], -1));
	printf("せん断弾性係数: %lf GPa\n", pow(In_plane_stiffness_matrix[2][2][99], -1) / height_sum);


	GyakuGyoretsu(Out_of_plane_stiffness_matrix, 99);
	printf("面外弾性特性: \n");
	printf("x方向のヤング率: %lf GPa\n", pow(Out_of_plane_stiffness_matrix[0][0][99], -1) * pow(height_sum, -3) * 12);
	printf("y方向のヤング率: %lf GPa\n", pow(Out_of_plane_stiffness_matrix[1][1][99], -1) * pow(height_sum, -3) * 12);
	printf("x方向のポアソン比: %lf\n", -Out_of_plane_stiffness_matrix[1][0][99] * pow(Out_of_plane_stiffness_matrix[0][0][99], -1));
	printf("y方向のポアソン比: %lf\n", -Out_of_plane_stiffness_matrix[0][1][99] * pow(Out_of_plane_stiffness_matrix[1][1][99], -1));
	printf("せん断弾性係数: %lf GPa\n", pow(Out_of_plane_stiffness_matrix[2][2][99], -1) * pow(height_sum, -3) * 12);


	/*double A[3][3], B[3][3], C[3][3];
	Input(A);
	Display(A);
	Input(B);
	Display(B);

	Gyoretsu_Kakesan(A, B, C);
	Display(C);
	GyakuGyoretsu(A);
	Display(A);*/

}
void Display(double Q[3][3][100], int I)//行列の表示
{
	// 表示
	int i, j;
	for (i = 0; i < 3; i++)
	{
		for (j = 0; j < 3; j++)
		{
			printf("%lf ", Q[i][j][I]);
		}
		printf("\n");
	}
}

int check(double mat[N][N][100], double inv[N][N][100], int I)
{

	double inner_product;
	int i, j, k;
	double ans;
	double diff;
	for (i = 0; i < N; i++)
	{
		for (j = 0; j < N; j++)
		{
			inner_product = 0;
			for (k = 0; k < N; k++)
			{
				inner_product += mat[i][k][I] * inv[k][j][I];
			}

			ans = (i == j) ? 1 : 0;
			diff = ans - inner_product;
			if (fabs(diff) > MAX_ERR) {
				return 0;
			}
		}
	}

	return 1;
}

double Tenchi(double Q[3][3])
{
	int i, j;
	double Q_a[3][3];
	for (i = 0; i < 3; i++)
	{
		for (j = 0; j < 3; j++)
		{
			Q_a[i][j] = Q[i][j];
		}
	}
	for (i = 0; i < 3; i++)
	{
		for (j = 0; j < 3; j++)
		{
			Q[j][i] = Q_a[i][j];
		}
	}
	return Q[3][3];
}

void Syokika(double Q[3][3][100], int I)
{
	int i, j;
	for (i = 0; i < 3; i++)
	{
		for (j = 0; j < 3; j++)
		{
			Q[i][j][I] = 0.0;
		}
	}
}

void GyakuGyoretsu(double Q[3][3][100], int I)
{
	/* 逆行列を求める行列用の２次元配列 */
	double mat[N][N][100];

	/* 逆行列用の２次元配列 */
	double inv[N][N][100];

	/* 掃き出し法を行う行列 */
	double sweep[N][N * 2][100];

	int i; /* 行 */
	int j; /* 列 */
	int k; /* 注目対角成分が存在する列 */

	double a; /* 定数倍用 */

	/* 正方行列の各成分をセット */
	for (i = 0; i < 3; i++)
	{
		for (j = 0; j < 3; j++)
		{
			mat[i][j][I] = Q[i][j][I];
		}
	}

	for (i = 0; i < N; i++)
	{
		for (j = 0; j < N; j++)
		{
			/* sweepの左側に逆行列を求める行列をセット */
			sweep[i][j][I] = mat[i][j][I];

			/* sweepの右側に単位行列をセット */
			sweep[i][N + j][I] = (i == j) ? 1 : 0;
		}
	}


	/* 全ての列の対角成分に対する繰り返し */
	for (k = 0; k < N; k++)
	{

		/* 最大の絶対値を注目対角成分の絶対値と仮定 */
		double max = fabs(sweep[k][k][I]);
		int max_i = k;

		/* k列目が最大の絶対値となる行を探す */
		for (i = k + 1; i < N; i++)
		{
			if (fabs(sweep[i][k][I]) > max)
			{
				max = fabs(sweep[i][k][I]);
				max_i = i;
			}
		}

		if (fabs(sweep[max_i][k][I]) <= MAX_ERR)
		{
			/* 逆行列は求められない */
			printf("逆行列は求められません...\n");
			return;
		}

		/* 操作（１）：k行目とmax_i行目を入れ替える */
		if (k != max_i) {
			for (j = 0; j < N * 2; j++)
			{
				double tmp = sweep[max_i][j][I];
				sweep[max_i][j][I] = sweep[k][j][I];
				sweep[k][j][I] = tmp;
			}
		}

		/* sweep[k][k]に掛けると1になる値を求める */
		a = 1 / sweep[k][k][I];

		/* 操作（２）：k行目をa倍する */
		for (j = 0; j < N * 2; j++)
		{
			/* これによりsweep[k][k]が1になる */
			sweep[k][j][I] *= a;
		}

		/* 操作（３）によりk行目以外の行のk列目を0にする */
		for (i = 0; i < N; i++)
		{
			if (i == k) {
				/* k行目はそのまま */
				continue;
			}

			/* k行目に掛ける値を求める */
			a = -sweep[i][k][I];

			for (j = 0; j < N * 2; j++)
			{
				/* i行目にk行目をa倍した行を足す */
				/* これによりsweep[i][k]が0になる */
				sweep[i][j][I] += sweep[k][j][I] * a;
			}
		}
	}

	/* sweepの右半分がmatの逆行列 */
	for (i = 0; i < N; i++)
	{
		for (j = 0; j < N; j++)
		{
			inv[i][j][I] = sweep[i][N + j][I];
		}
	}

	/*Qにinvを代入*/
	for (i = 0; i < N; i++)
	{
		for (j = 0; j < N; j++)
		{
			Q[i][j][I] = inv[i][j][I];
		}
	}

	/* 逆行列invを表示
	for (i = 0; i < N; i++)
	{
		for (j = 0; j < N; j++)
		{
			printf("%f, ",Q[i][j]);
		}
		printf("\n");
	}*/

	/* 検算
	if (check(mat, inv))
	{
		printf("invはmatの逆行列です！！\n");
	}
	else
	{
		printf("invはmatの逆行列になってません...\n");
	}*/
}

void Gyoretsu_Kakesan(double Q[3][3][100], double R[3][3][100], double result[3][3][100], int I)/*この行列の積（演算結果）は最後に用意された行列に入る事に注意*/
{
	int i, j, k;
	Syokika(result, I);

	for (i = 0; i < 3; i++)
	{
		for (j = 0; j < 3; j++)
		{
			for (k = 0; k < 3; k++)
			{

				result[i][j][I] += Q[i][k][I] * R[k][j][I];

			}
		}
	}
}

void Input(double A[3][3][100], int I)
{

	int i, j;
	Syokika(A, I);
	for (i = 0; i < 3; i++)
	{
		for (j = 0; j < 3; j++)
		{
			printf("A_%d%dを入力してください。: ", i, j);
			scanf_s("%lf", &A[i][j][I]);
		}
	}
}

void Stress_Strain(double T_sigma[3][3][100], double Q[3][3][100], double T_epsilon[3][3][100], double result[3][3][100], int I)
{
	//printf("x-y座標系の応力とひずみの関係式の係数行列\n=====================\n");
	GyakuGyoretsu(T_sigma, I);
	double result1[3][3][100];
	Syokika(result1, I);
	Gyoretsu_Kakesan(T_sigma, Q, result1, I);
	Gyoretsu_Kakesan(result1, T_epsilon, result, I);
	//Display(result);
	//printf("=====================\n");
}


//下三つでの操作でIに代入可能なのは1からn-1までであることに注意する。
void In_plane_stiffness_matrix_function(double Q_bar[3][3][100], double Thickness[100], double In_plane_stiffness_matrix[3][3][100], int I)
{
	/*面内剛性行列(In_plane_stiffness_matrix)を求める*/

		int i, j;
		{
			for (i = 0; i < 3; i++)
			{
				for (j = 0; j < 3; j++)
				{
					In_plane_stiffness_matrix[i][j][99] += Q_bar[i][j][I] * (Thickness[I + 1] - Thickness[I]);
				}
			}
		}

}

void Tensile_bending_coupling_matrix_function(double Q_bar[3][3][100], double Thickness[100], double Tensile_bending_coupling_matrix[3][3][100], int I)
{
	/*引っ張り-曲げカップリング行列を求める*/
	int i, j;

	{
		for (i = 0; i < 3; i++)
		{
			for (j = 0; j < 3; j++)
			{
				Tensile_bending_coupling_matrix[i][j][99] += Q_bar[i][j][I] * (pow(Thickness[I + 1], 2) - pow(Thickness[I], 2)) * 0.5;
			}
		}
	}
}

void Out_of_plane_stiffness_matrix_function(double Q_bar[3][3][100], double Thickness[100], double Out_of_plane_stiffness_matrix[3][3][100], int I)
{
	/*面外剛性行列を求める*/
	int i, j;
	{
		for (i = 0; i < 3; i++)
		{
			for (j = 0; j < 3; j++)
			{
				Out_of_plane_stiffness_matrix[i][j][99] += Q_bar[i][j][I] * (pow(Thickness[I + 1], 3) - pow(Thickness[I], 3)) * 0.333333;
			}
		}
	}
}
// プログラムの実行: Ctrl + F5 または [デバッグ] > [デバッグなしで開始] メニュー
// プログラムのデバッグ: F5 または [デバッグ] > [デバッグの開始] メニュー

// 作業を開始するためのヒント: 
//    1. ソリューション エクスプローラー ウィンドウを使用してファイルを追加/管理します 
//   2. チーム エクスプローラー ウィンドウを使用してソース管理に接続します
//   3. 出力ウィンドウを使用して、ビルド出力とその他のメッセージを表示します
//   4. エラー一覧ウィンドウを使用してエラーを表示します
//   5. [プロジェクト] > [新しい項目の追加] と移動して新しいコード ファイルを作成するか、[プロジェクト] > [既存の項目の追加] と移動して既存のコード ファイルをプロジェクトに追加します
//   6. 後ほどこのプロジェクトを再び開く場合、[ファイル] > [開く] > [プロジェクト] と移動して .sln ファイルを選択します
