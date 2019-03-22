#include<iostream>
#include<stdio.h>
#include <sstream>
#include <string>
#include <vector>
#include <stdlib.h>
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

void CreatSimplexTableB(vector < vector<float> > &SimplexTable, vector<float> &Z, vector<float> &CN, vector<float> &CB, vector<string> &bv, vector<string> &nbv, int n, int m, int slack_var, vector < vector<float> > &NewMatrix, float *b, float *Optim, float d)
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


void CreatSimplexTableA(vector < vector<float> > &SimplexTable, vector<float> &Z, vector<float> &CN, vector<float> &CB, vector<string> &bv, vector<string> &nbv, int n, int m, int slack_var, vector < vector<float> > &NewMatrix, float *b, float *Optim, float d)
{
	int i,j;
	float sum;
	for (i=0;i<n;i++)
	{
		ostringstream ss;
		ss << i+1+m;
		bv[i] = "x" + ss.str();
	}
	bv[n] = "Z ";
		
	for (i=0;i<m;i++)
	{
		ostringstream ss;
		ss << i+1;
		nbv[i] = "x" + ss.str();
	}
	
	for (i=0;i<m;i++)
	{
		for (j=0;j<n;j++)
		{
			SimplexTable[j][i] = -NewMatrix[i][j];
		}
	}
	for (i=0;i<n;i++)
	{
		SimplexTable[i][m] = -1;
	}
	
	for (i=0;i<m;i++)
	{
		CN[i] = -1;
	}
	
	for (i=0;i<n;i++)
	{
		CB[i] = 0;
	}
	
	for (i=0;i<m;i++)
	{	
		sum = 0;	
		for (j=0;j<n;j++)
		{
			sum = sum + CB[j] * SimplexTable[j][i];
		}
		sum = sum - CN[i];
		Z[i] = sum;
	}
	sum = 0;
	for (i=0;i<m;i++)
	{
		sum = sum + CB[i] * (-1);
	}
	Z[m] = sum;
}

void PrimalSimplex(int n, int m, vector < vector<float> > &STable, vector<string> &bv, vector<string> &nbv, vector<float> &Z, vector<float> &CN, vector<float> &CB, bool &condn, vector<float> &ratio)
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

void DualSimplex(int n, int m, vector < vector<float> > &STable, vector<string> &bv, vector<string> &nbv, vector<float> &Z, vector<float> &CN, vector<float> &CB, bool &condn, vector<float> &ratio)
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
		if (Z[i] != 0 and STable[k][i] != 0) ratio[i] = Z[i]/STable[k][i];
		else ratio[i] = 1/0.0000000001;
	}
	
	for (i=0;i<n;i++)
	{
		for (j=i;j<n;j++)
		{
			if (ratio[i] < ratio[p] & ratio[i] != 0 & STable[k][i] < 0)
			{
				p = i;
			}
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

void print_simplex(int n, int m, vector < vector<float> > &STable, vector<string> &bv, vector<string> &nbv, vector<float> &Z, vector<float> &CN, vector<float> &CB)
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

void Stablity(vector < vector<float> > &NewMatrix, int n, int m, bool &condn)
{
	int i,j,k=0,l=0;
	float r,c,minA[m],maxB[n];
	
	for (i=0;i<m;i++)
	{
		minA[i] = NewMatrix[i][0];
	}
	for (i=0;i<n;i++)
	{
		maxB[i] = NewMatrix[0][i];
	}
	
	for (i=0;i<m;i++)
	{
		for (j=0;j<n;j++)
		{
			if (NewMatrix[i][j] < minA[i]) minA[i] = NewMatrix[i][j];
		}
	}
	for (j=0;j<n;j++)
	{
		for (i=0;i<m;i++)
		{
			if (NewMatrix[i][j] > maxB[j]) maxB[j] = NewMatrix[i][j];
		}
	}
	
	r = minA[0];
	c = maxB[0];
	for (i=0;i<m;i++)
	{
		if (minA[i] > r) {r = minA[i]; k=i;}
	}
	for (i=0;i<n;i++)
	{
		if (maxB[i] < c) {c = maxB[i]; l=i;}
	}
	
	if (r == c) {
		cout << "Game is Stable" << endl; 
		condn = true;
		cout << "Game Value is " << r << endl;
		cout << "Strategy of A is " << k << endl;
		cout << "Strategy of B is " << l << endl;}
	else cout << "Game is Unstable" << endl;
}

void SolveDual(int n, int m, vector < vector<float> > &STable, vector<string> &bv, vector<string> &nbv, vector<float> &Z, vector<float> &CN, vector<float> &CB, bool &condn, vector<float> &ratio)
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

void SolvePrimal(int n, int m, vector < vector<float> > &STable, vector<string> &bv, vector<string> &nbv, vector<float> &Z, vector<float> &CN, vector<float> &CB, bool &condn, vector<float> &ratio)
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

void SolutionB(float SolB[], int n, float ProbB[],float &game_valueB, float min_value)
{
	int i,j,k;
	float sum=0;
	for (i=0;i<n;i++)
	{
		sum = sum + SolB[i];
	}
	game_valueB = 1/sum;
	for (i=0;i<n;i++)
	{
		ProbB[i] = SolB[i] * game_valueB;
	}
	game_valueB = game_valueB + min_value;
	
	cout << "Probabilities of Strategies used by B are : ";
	for (i=0;i<n;i++)
	{
		cout << ProbB[i] << "    ";
	}
	cout << endl;
}

void SolutionA(float SolA[], int m, float ProbA[],float &game_valueA, float min_value)
{
	int i,j,k;
	float sum=0;
	for (i=0;i<m;i++)
	{
		sum = sum + SolA[i];
	}
	game_valueA = 1/sum;
	for (i=0;i<m;i++)
	{
		ProbA[i] = SolA[i] * game_valueA;
	}
	game_valueA = game_valueA + min_value;
	
	cout << "Probabilities of Strategies used by A are : ";
	for (i=0;i<m;i++)
	{
		cout << ProbA[i] << "    ";
	}
	cout << endl;
}

int main()
{
	int n, m, max,ith=0,NI;
  	int i,j,k,slack_var=0;
  	bool condn = false;
  	cout << "Enter values of n (Strategies of B) and m(Strategies of A) " << endl;
  	cin >> n >> m;

	// For inequalities system
	int a[m];
	float min_value;
	for (i=0;i<m;i++)
	{
		a[i] = 1;
	}
	
  	float A[m][n], b[m], Optim[n],d = 0;
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
  	
  	min_value = 0;
  	for (i=0;i<m;i++)
  	{
  		for (j=0;j<n;j++)
  		{
  			if (min_value > A[i][j]) min_value = A[i][j];
  		}
  	}
  	
  	if (min_value < 0)
  	{
  		for (i=0;i<m;i++)
  		{
  			for (j=0;j<n;j++)
  			{
  				A[i][j] = A[i][j] - min_value;
  			}
  		}
  	}

  	for (i=0;i<n;i++)
  	{
  		Optim[i] = 1;
  	}
  	
  	for (i=0;i<m;i++)
  	{
  		b[i] = 1;
  	}

	int ns = n+slack_var,g;
	vector < vector<float> > NewMatrix(m, vector<float> (ns, 0));

	for (i=0;i<m;i++)
	{
		for (j=0;j<n;j++)
		{
			NewMatrix[i][j] = A[i][j];
		}
	}
	
	Stablity(NewMatrix,n,m,condn);
	
	vector < vector<float> > SimplexTable(m, vector<float> (n+1, 0));
	vector<float> Z(n+1), CN(n), CB(m), ratio(n);
	
	vector<string> bv(slack_var+1), nbv(n);
	
	cout << "Table" << endl; 
	for (i=0;i<m;i++)
	{	
		for (j=0;j<n;j++)
		{
			cout << NewMatrix[i][j] << "    ";
		}
		cout << endl;
	}
	
	if (!condn)
	{
		// Solving for A
		SimplexTable.resize(n+m,vector<float>(m+n+1));
		Z.resize(m+1);
		CN.resize(m);
		CB.resize(n);
		ratio.resize(m);
		bv.resize(n+1);
		nbv.resize(m);
		CreatSimplexTableA(SimplexTable, Z, CN, CB, bv, nbv, n, m, slack_var, NewMatrix, b, Optim, d);
	
		cout << "Solving for A:" << endl;
		print_simplex(m, n, SimplexTable, bv, nbv, Z, CN, CB);
		cout << endl;
		
		SolveDual(m, n, SimplexTable, bv, nbv, Z, CN, CB, condn, ratio);
		
		SolvePrimal(m, n, SimplexTable, bv, nbv, Z, CN, CB, condn, ratio);
		cout << "Solution table for A: " << endl;
		print_simplex(m, n, SimplexTable, bv, nbv, Z, CN, CB);
		cout << endl;
		
		float SolA[m],ProbA[m],game_valueA;
		int idx=0;
		string temp_s;
		for (i=0;i<m;i++) SolA[i] = 0;
		for (i=0;i<n;i++)
		{
			temp_s = bv[i][1];
			idx = atoi(temp_s.c_str());
			if (idx < (n+1)) SolA[idx-1] = SimplexTable[i][m];
		}
		
		SolutionA(SolA,m,ProbA, game_valueA,min_value);
		
		cout << "Game Value (A) : " << game_valueA << endl;
		
		// Solving for B
		SimplexTable.resize(m,vector<float>(n+1));
		Z.resize(n+1);
		CN.resize(n);
		CB.resize(m);
		ratio.resize(n);
		bv.resize(m+1);
		nbv.resize(n);
		CreatSimplexTableB(SimplexTable, Z, CN, CB, bv, nbv, n, m, slack_var, NewMatrix, b, Optim, d);
		cout << "Solving for Player B :" << endl;
		print_simplex(n, m, SimplexTable, bv, nbv, Z, CN, CB);
		cout << endl;
		
		SolvePrimal(n, m, SimplexTable, bv, nbv, Z, CN, CB, condn, ratio);
		cout << "Solution for B: " << endl;	
		print_simplex(n, m, SimplexTable, bv, nbv, Z, CN, CB);
		cout << endl;
		
		float SolB[n],ProbB[n],game_valueB;
		//int idx=0;
		//string temp_s;
		for (i=0;i<n;i++) SolB[i] = 0;
		for (i=0;i<m;i++)
		{
			temp_s = bv[i][1];
			idx = atoi(temp_s.c_str());
			if (idx < (n+1)) SolB[idx-1] = SimplexTable[i][n];
		}
		
		SolutionB(SolB,n,ProbB, game_valueB,min_value);
		
		cout << "Game Value (B) : " << game_valueB << endl;
	}
		
  	return 0;
}
