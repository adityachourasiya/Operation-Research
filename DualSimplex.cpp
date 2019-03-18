#include<iostream>
#include<stdio.h>
#include <sstream>
#include <string>
#include <vector>
using namespace std;

int cbps = 0;
int number_Iteration = 0;

float round(float var)
{
    	float value = (int)(var * 100 + .5);
    	return (float)value / 100;
} 

void combinationUtil(vector < vector<int> > &comb, int arr[], int data[], int start, int end, int index, int r);

void printCombination(vector < vector<int> > &comb, int arr[], int n, int r) 
{ 
    	int data[r]; 
    	combinationUtil(comb,arr, data, 0, n-1, 0, r); 
} 

void p(vector < vector<int> > &comb, int data[], int r)
{
	int i;
	for (i=0;i<r;i++)
	{			
		comb[cbps][i] = data[i];
	}
	cbps++;
}

void combinationUtil(vector < vector<int> > &comb, int arr[], int data[], int start, int end, int index, int r) 
{  
    	if (index == r) 
    	{ 
        	p(comb,data,r);
        	return; 
    	} 

    	for (int i=start; i<=end && end-i+1 >= r-index; i++) 
    	{ 
        	data[index] = arr[i];
        	combinationUtil(comb, arr, data, i+1, end, index+1, r); 
    	} 
} 

void gaussian(int n, int m, vector < vector<float> > &A, float *b, int *newtemp, int ns, vector < vector<float> > &basic_sols,int len)
{
	// Gaussian Elemination Method
  	int i,j,k;
  	
  	for (int k=0;k<n-1;k++)
  	{
		b[k] /= A[k][k];
		for (i=n-1;i>=0;i--)
		{			
			A[k][i] /= A[k][k];			
		}	
  		for (i=k+1;i<m;i++)
  		{  
			b[i] = b[i] - (A[i][k] * b[k]);		
  			for (j=n-1;j>=0;j--)
  			{
  				A[i][j] = A[i][j] - (A[i][k] * A[k][j]);
  			}
  		}		
  	}
	
	if (A[m-1][n-1] != 0)
	{
		b[m-1] = b[m-1]/(A[m-1][n-1]);
		A[m-1][n-1] /= (A[m-1][n-1]);
	}


	for (i=0;i<m-1;i++)
	{
		for (j=i+1;j<m;j++)
		{
			b[i] = b[i] - (A[i][j] * b[j]);
			for (k=n-1;k>=0;k--)
			{
				A[i][k] = A[i][k] - (A[i][j] * A[j][k]);				
			}
		}
	}	
	
	j=0;
	for (i=0;i<ns;i++)
	{
		if (newtemp[j] == i)
		{
			basic_sols[len][i] = b[j];
			j++;
		}
		else
		{
			basic_sols[len][i] =  0;
		}
	}	
  		
}

int fact(int n){
    	if(n==0) return 1;
    	if (n>0) return n*fact(n-1);
}

int NCR(int n,int r){
    	if(n==r) return 1;
    	if (r==0&&n!=0) return 1;
    	else return (n*fact(n-1))/(fact(r)*fact(n-r));
}

void CreatSimplexTable(vector < vector<float> > &SimplexTable, float *Z, float *CN, float *CB, vector<string> &bv, vector<string> &nbv, int n, int m, int slack_var, vector < vector<float> > &NewMatrix, float *b, float *Optim, float d)
{
	int i,j;
	float sum;
	for (i=0;i<slack_var;i++)
	{
		ostringstream ss;
		ss << i+1+n;
		bv[i] = "x" + ss.str();
	}
		
	for (i=0;i<n;i++)
	{
		ostringstream ss;
		ss << i+1;
		nbv[i] = "x" + ss.str();
	}
	bv[slack_var] = "Z ";
	
	for (i=0;i<m;i++)
	{
		for (j=0;j<n;j++)
		{
			SimplexTable[i][j] = NewMatrix[i][j];
		}
		SimplexTable[i][j] = b[i];
	}
	
	for (i=0;i<n;i++)
	{
		CN[i] = Optim[i];
	}
	
	for (i=0;i<m;i++)
	{
		CB[i] = 0;
	}
	
	for (i=0;i<n;i++)
	{	
		sum = 0;	
		for (j=0;j<m;j++)
		{
			sum = sum + CB[j] * SimplexTable[j][i];
		}
		sum = sum - CN[i];
		Z[i] = sum;
	}
	sum = 0;
	for (i=0;i<m;i++)
	{
		sum = sum + CB[i] * b[i];
	}
	Z[n] = sum;
}

void PrimalSimplex(int n, int m, vector < vector<float> > &STable, vector<string> &bv, vector<string> &nbv, float Z[], float CN[], float CB[], bool &condn, float *ratio)
{
	number_Iteration++;
	int i, j, k=0,p=0;
	float pivot, sum, x;
	string temp;
	for (i=0;i<n;i++)
	{
		for (j=i;j<n;j++)
		{
			if (Z[i] < Z[k])
			{
				k = i;
			}
		}
	}
	
	for (i=0;i<m;i++)
	{
		if (STable[i][k] > 0) ratio[i] = STable[i][n]/STable[i][k];
		else ratio[i] = 0;
	}
	
	for (i=0;i<m;i++)
	{
		for (j=i;j<m;j++)
		{
			if (ratio[i] > ratio[p])
			{
				p = i;
			}
		}
	}
	
	for (i=0;i<m;i++)
	{
		for (j=i;j<m;j++)
		{
			if (ratio[i] < ratio[p] & ratio[i] > 0)
			{
				p = i;
			}
		}
	}
	
	pivot = STable[p][k];
	
	for (i=0;i<m;i++)
	{
		if (i != p)
		{
			for (j=0;j<n+1;j++)
			{
				if (j != k)
				{
					STable[i][j] = ((STable[i][j] * STable[p][k]) - (STable[i][k] * STable[p][j])) / pivot;
				}
			}
		}
	}
	
	for (i=0;i<n+1;i++)
	{
		if (i != k)
		{
			STable[p][i] = STable[p][i]/pivot;
		}
	}
	for (i=0;i<m;i++)
	{
		if (i != p)
		{
			STable[i][k] = -(STable[i][k]/pivot);
		}
	}	
	STable[p][k] = 1/STable[p][k];
	
	x = CN[k];
	CN[k] = CB[p];
	CB[p] = x;
	
	temp = bv[p];
	bv[p] = nbv[k];
	nbv[k] = temp;
	
	for (i=0;i<n;i++)
	{
		sum = 0;
		for (j=0;j<m;j++)
		{
			sum += (STable[j][i] * CB[j]);
		}
		Z[i] = sum - CN[i];
	}
	sum = 0;
	for (j=0;j<m;j++)
	{
		sum += (STable[j][i] * CB[j]);
	}
	Z[n] = sum;	
	
	for (i=0;i<m;i++)
	{
		if (STable[i][n] < 0)
		{
			condn = true;
		}
	}	
}

void DualSimplex(int n, int m, vector < vector<float> > &STable, vector<string> &bv, vector<string> &nbv, float Z[], float CN[], float CB[], bool &condn, float *ratio)
{
	number_Iteration++;
	int i, j, k=0,p=0;
	float pivot, sum, x;
	string temp;
	for (i=0;i<m;i++)
	{
		for (j=i;j<m;j++)
		{
			if (STable[i][n] < STable[k][n])
			{
				k = i;
			}
		}
	}
	
	for (i=0;i<n;i++)
	{
		if (Z[i] != 0 and STable[k][i] != 0) ratio[i] = -(Z[i]/STable[k][i]);
		else ratio[i] = 1/0.0000000001;
	}
	
	for (i=0;i<n;i++)
	{
		if (ratio[i] < ratio[p] & ratio[i] != 0 & STable[k][i] < 0)
		{
			p = i;
		}
	}
	
	pivot = STable[k][p];
	
	for (i=0;i<m;i++)
	{
		if (i != k)
		{
			for (j=0;j<n+1;j++)
			{
				if (j != p)
				{
					STable[i][j] = ((STable[i][j] * STable[k][p]) - (STable[i][p] * STable[k][j])) / pivot;
				}
			}
		}
	}
	
	for (i=0;i<n+1;i++)
	{
		if (i != p)
		{
			STable[k][i] = STable[k][i]/pivot;
		}
	}
	for (i=0;i<m;i++)
	{
		if (i != k)
		{
			STable[i][p] = -(STable[i][p]/pivot);
		}
	}	
	STable[k][p] = 1/STable[k][p];
	
	x = CN[p];
	CN[p] = CB[k];
	CB[k] = x;
	
	temp = bv[k];
	bv[k] = nbv[p];
	nbv[p] = temp;
	
	for (i=0;i<n;i++)
	{
		sum = 0;
		for (j=0;j<m;j++)
		{
			sum += (STable[j][i] * CB[j]);
		}
		Z[i] = sum - CN[i];
	}
	sum = 0;
	for (j=0;j<m;j++)
	{
		sum += (STable[j][i] * CB[j]);
	}
	Z[n] = sum;		
}

void IntegerProg(int n, int m, vector < vector<float> > &STable, vector<string> &bv, vector<string> &nbv, float Z[], float CN[], float CB[], bool &condn, float *ratio)
{
	int i,j,k=0;
	float fraction[m],frac,sum;
	for (i=0;i<m;i++)
	{
		if (STable[i][n] > 0)
		{
			frac = STable[i][n] - (int)(STable[i][n]+0.01);
			fraction[i] = frac;
		}
		else
		{
			frac = STable[i][n] - (int)(STable[i][n]+0.01) + 1;
			fraction[i] = frac;
		}
		
	}
	
	for (i=0;i<m;i++)
	{
		if (fraction[k] < fraction[i]) k = i;
	}
	
	for (i=0;i<n;i++)
	{
		if (STable[k][i] > 0)
		{
			STable[m][i] = -(STable[k][i] - (int)STable[k][i]);
		}
		else
		{
			STable[m][i] = -(STable[k][i] - (int)STable[k][i] + 1);
		}
	}
	
	m++;
	STable[m-1][n] = -(STable[k][n] - (int)STable[k][n]);
	bv[m] = bv[m-1];
	ostringstream ss;
	ss << m+n;
	bv[m-1] = "x" + ss.str();
	for (i=0;i<n;i++)
	{
		sum = 0;
		for (j=0;j<m;j++)
		{
			sum += (STable[j][i] * CB[j]);
		}
		Z[i] = sum - CN[i];
	}
	sum = 0;
	for (j=0;j<m;j++)
	{
		sum += (STable[j][i] * CB[j]);
	}
	Z[n] = sum;
}

void print_simplex(int n, int m, vector < vector<float> > &STable, vector<string> &bv, vector<string> &nbv, float Z[], float CN[], float CB[])
{
	int i,j,k;
	cout << "        CN   |";
	for (i=0;i<n;i++)
	{
		if (CN[i] < 0) printf("     %.2f |", round(CN[i]));
		else printf("      %.2f |", round(CN[i]));
	}
	cout << endl;
	cout << "CB  |";
	cout << "  BV\\NBV|";
	for (i=0;i<n;i++)
	{
		printf("       %s  |", nbv[i].c_str());
	}
	cout << "     Solution|" << endl;
	
	
	for (i=0;i<m;i++)
	{
		printf("%.2f|", round(CB[i]));
		printf("     %s |", bv[i].c_str());
		for (j=0;j<n+1;j++)
		{
			if (STable[i][j] < 0) printf("     %.2f |", round(STable[i][j]));
			else printf("      %.2f |", round(STable[i][j]));
		}
		cout << endl;
	}
	
	cout << "           Z |";
	for (i=0;i<n+1;i++)
	{
		if (Z[i] < 0) printf("     %.2f |", round(Z[i]));
		else printf("      %.2f |", round(Z[i]));
	}
	cout << endl;
}

void SolveDual(int n, int m, vector < vector<float> > &STable, vector<string> &bv, vector<string> &nbv, float Z[], float CN[], float CB[], bool &condn, float *ratio)
	{
		int i,j,k;
		for (i=0;i<m;i++)
		{
			if (STable[i][n] < 0)
			{	
				DualSimplex(n, m, STable, bv, nbv, Z, CN, CB, condn, ratio);
				i=-1;
			}
			if (condn == true)
			{
				break;
			}
		}
	}

void SolvePrimal(int n, int m, vector < vector<float> > &STable, vector<string> &bv, vector<string> &nbv, float Z[], float CN[], float CB[], bool &condn, float *ratio)
	{
		int i,j,k;
		for (i=0;i<n;i++)
		{
			if (Z[i] < 0)
			{
				PrimalSimplex(n, m, STable, bv, nbv, Z, CN, CB, condn, ratio);
				i=-1;
			}
			if (condn == true)
			{
				break;
			}
		}
	}

int main()
{
	int n, m, max,ith=0,NI;
  	int i,j,k,slack_var=0;
  	bool condn = false;
  	cout << "Enter values of n (number of unknowns) and m(number of equations) " << endl;
  	cin >> n >> m;

	// For inequalities system
	int a[m];
	cout << "Enter the sign of equations 0 (if =), -1 (if >=) or 1 (if <=)" << endl;
	for (i=0;i<m;i++)
	{
		cin >> a[i];
	}
	
	
  	float A[m][n], b[m], Optim[n],d;
  	// Input A matrix
  	cout << "Enter the elements rowise" << endl;
  	for (i=0;i<m;i++)
  	{
  		for (j=0;j<n;j++)
  		{
  			cin >> A[i][j];
  		}
		slack_var++;
  	}

  	// Inoput b matrix
  	cout << "Enter the b matrix (shape of b matrix must be m X 1)" << endl;
  	for (i=0;i<m;i++)
  	{
  		cin >> b[i];
  	}
  	
  	cout << "Define your optimization problem: Enter 0 (Minimizaion) or 1 (Maximization)" << endl;
  	cin >> max;
  	cout << "Enter the coeff of your optimization problem in order: (length should be equal to number of unknown given) add 0 if there is no constant at the end" << endl;
  	for (i=0;i<n;i++)
  	{
  		cin >> Optim[i];
  	}
  	cin >> d;	
  	
  	for (i=0;i<n;i++)
  	{
  		if (max == 0)
  		{
  			Optim[i] = Optim[i]*(-1);
  		}
  	}
  	
  	for (i=0;i<m;i++)
  	{
  		if (a[i] == -1)
  		{
  			b[i] = b[i]*(-1);
  		}
  	}

	int ns = n+slack_var,g;
	vector < vector<float> > NewMatrix(m, vector<float> (ns, 0));
	cout << "Enter \n 1 List of all BFS \n 2 Number of Iteration \n 3 List of all Non-basic variables along with net evaluations in ith iteration \n 4 List of Basic variables along with min ratios in ith iteration \n 5 simplex table of ith iteration \n 6 Optimal solution \n 7 Exit" << endl;
	cin >> g;

	// Making NewMatrix = A + slack variables

	for (i=0;i<m;i++)
	{
		if (a[i] == -1)
		{
			for (j=0;j<n;j++)
			{
				NewMatrix[i][j] = -A[i][j];
			}
		}
		else
		{
			for (j=0;j<n;j++)
			{
				NewMatrix[i][j] = A[i][j];
			}
		}
	}
	
	for (i=0;i<m;i++)
	{
		if (a[i] == 0)
		{
			for (j=n;j<ns;j++)
			{
				NewMatrix[i][j] = 0;
			}			
		}
		else if (a[i] == 1)
		{
			for (j=0;j<slack_var;j++)
			{
				if (i == j)
				{
					NewMatrix[i][j+n] = a[i];
				}
				else
				{
					NewMatrix[i][j+n] = 0;
				}
			}
		}

		else if (a[i] == -1)
		{
			for (j=0;j<slack_var;j++)
			{
				if (j == i)
				{
					NewMatrix[i][j+n] = a[i];
				}
				else
				{
					NewMatrix[i][j+n] = 0;
				}
			}
		}
		else
		{
			cout << "Error in providing the input" << endl;
		}
	}
	
	vector < vector<float> > SimplexTable(m, vector<float> (n+1, 0));
	float Z[n+1], CN[n], CB[m], ratio[n];
	
	vector<string> bv(slack_var+1), nbv(n);
	
	CreatSimplexTable(SimplexTable, Z, CN, CB, bv, nbv, n, m, slack_var, NewMatrix, b, Optim, d);
	print_simplex(n, m, SimplexTable, bv, nbv, Z, CN, CB);
	cout << endl;
	
	SolveDual(n, m, SimplexTable, bv, nbv, Z, CN, CB, condn, ratio);
	
	print_simplex(n, m, SimplexTable, bv, nbv, Z, CN, CB);
	cout << endl;
	
	SolvePrimal(n, m, SimplexTable, bv, nbv, Z, CN, CB, condn, ratio);
	
	print_simplex(n, m, SimplexTable, bv, nbv, Z, CN, CB);
	cout << "starts" << endl;
	
	float frac;
	for (i=0;i<m;i++)
	{
		frac = SimplexTable[i][n] - (int)(SimplexTable[i][n]);
		if (frac > 0 & frac < 1)
		{
			SimplexTable.resize(m+1,vector<float>(n+1));
			bv.resize(m+2);
			IntegerProg(n, m, SimplexTable, bv, nbv, Z, CN, CB, condn, ratio);
			m++;
			SolveDual(n, m, SimplexTable, bv, nbv, Z, CN, CB, condn, ratio);
			i=-1;
		}
	}
	
	cout << "ends" << endl;
	cout << "1" << endl;
	SimplexTable.resize(m+1,vector<float>(n+1));
	bv.resize(m+1+1);
	cout << "Resizing" << endl;
	print_simplex(n, m, SimplexTable, bv, nbv, Z, CN, CB);
	cout << endl;
	IntegerProg(n, m, SimplexTable, bv, nbv, Z, CN, CB, condn, ratio);
	m++;
	print_simplex(n, m, SimplexTable, bv, nbv, Z, CN, CB);
	cout << endl;
	cout << "Dual" << endl;
	SolveDual(n, m, SimplexTable, bv, nbv, Z, CN, CB, condn, ratio);
	print_simplex(n, m, SimplexTable, bv, nbv, Z, CN, CB);
	cout << endl;
	cout << "Primal" << endl;
	SolvePrimal(n, m, SimplexTable, bv, nbv, Z, CN, CB, condn, ratio);
	print_simplex(n, m, SimplexTable, bv, nbv, Z, CN, CB);
	cout << endl;
	
	for (i=0;i<m;i++)
	{
		frac = SimplexTable[i][n] - (int)(round(SimplexTable[i][n]));
		if (round(frac) != 0)
		{
			cout << round(frac) << endl;
		}
	}
	
	cout << "2" << endl;
	SimplexTable.resize(m+1,vector<float>(n+1));
	bv.resize(m+1+1);
	cout << "Resizing" << endl;
	print_simplex(n, m, SimplexTable, bv, nbv, Z, CN, CB);
	cout << endl;
	IntegerProg(n, m, SimplexTable, bv, nbv, Z, CN, CB, condn, ratio);
	m++;
	print_simplex(n, m, SimplexTable, bv, nbv, Z, CN, CB);
	cout << endl;
	cout << "Dual" << endl;
	SolveDual(n, m, SimplexTable, bv, nbv, Z, CN, CB, condn, ratio);
	print_simplex(n, m, SimplexTable, bv, nbv, Z, CN, CB);
	cout << endl;
	cout << "Primal" << endl;
	SolvePrimal(n, m, SimplexTable, bv, nbv, Z, CN, CB, condn, ratio);
	print_simplex(n, m, SimplexTable, bv, nbv, Z, CN, CB);
	cout << endl;
	
	for (i=0;i<m;i++)
	{
		frac = SimplexTable[i][n] - (int)(round(SimplexTable[i][n]));
		if (round(frac) != 0)
		{
			cout << -(int)(round(SimplexTable[i][n])) << endl;
			cout << i << "   " << round(frac) << endl;
		}
	}
	
		
  	return 0;
}
