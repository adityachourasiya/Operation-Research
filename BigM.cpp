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

void p(vector < vector <int> > &comb, int data[], int r)
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

void CreatSimplexTable(vector < vector<float> > &SimplexTable, float *Z, float *CN, float *CB, string *bv, string *nbv, int n, int m, int slack_var, vector < vector<float> > &NewMatrix, float *b, float *Optim, float d, int *a, float M)
{
	int i,j;
	float sum=0;
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
		if (a[i] == -1 | a[i] == 0)
		{
			CB[i] = -M;
			Z[n-i-1] = M;
		}
		else CB[i] = 0;
	}
	
	for (i=0;i<n;i++)
	{
		Z[i] = -Optim[i];
	}
	sum = 0;
	for (j=0;j<m;j++)
	{
		sum += (SimplexTable[j][n] * CB[j]);
	}
	Z[n] = sum+d;
}

void print_simplex(int n, int m, vector < vector<float> > &STable, string bv[], string nbv[], float Z[], float CN[], float CB[])
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


void Simplex(int n, int m, vector < vector<float> > &STable, string bv[], string nbv[], float Z[], float CN[], float CB[], bool &condn, float *ratio)
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



int main()
{
	int n, m, max,ith=0,NI,surplus_var = 0,g;
  	int i,j,k,slack_var=0;
  	bool condn = false, alternate = true;
  	cout << "Enter values of n (number of unknowns) and m(number of equations) " << endl;
  	cin >> n >> m;

	// For inequalities system
	int a[m],a_temp[m];
	cout << "Enter the sign of equations 0 (if =), -1 (if >=) or 1 (if <=)" << endl;
	for (i=0;i<m;i++)
	{
		cin >> a[i];
	}
	
	
  	float A[m][n], b[m], Optim[n],d, M;
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
  	
  	cout << "Enter big M value" << endl;
  	cin >> M;	

	cout << "Enter \n 1 List of all BFS \n 2 Number of Iteration \n 3 List of all Non-basic variables along with net evaluations in ith iteration \n 4 List of Basic variables along with min ratios in ith iteration \n 5 simplex table of ith iteration \n 6 Optimal solution \n 7 Exit" << endl;
	cin >> g;
	
	for (i=0;i<m;i++)
	{
		if (a[i] == -1)
		{
			surplus_var += 1;
		}
	}
		
	
	int ns = n+slack_var+surplus_var;
	vector < vector<float> > NewMatrix(m, vector<float> (ns, 0));
	float NewB[ns];

	// Making NewMatrix = A + slack variables

	for (i=0;i<m;i++)
	{
		for (j=0;j<n;j++)
		{
			NewMatrix[i][j] = A[i][j];
		}
		NewB[i] = b[i];
	}
	
	for (i=0;i<ns-m;i++)
	{
		NewB[i+m] = 0;
	}
	
	for (i=0;i<surplus_var;i++)
	{
		for (j=0;j<m;j++)
		{
			NewMatrix[j][n+i] = 0;
		}
	}
	
	j=n;
	for (i=0;i<m;i++)
	{
		if (a[i] == -1)
		{
			NewMatrix[i][j] = -1;
			j++;
		}
	}
	
	n += surplus_var;

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
	
	// Chosing m square matrix at a time

	int bps = NCR(ns,m),x;
	int arr[ns];
	vector < vector<int> > comb(bps, vector<int> (m, 0));

	for (i=0;i<ns;i++)
	{
		arr[i] = i;
	}

	for (i=0;i<bps;i++)
	{
		for (j=0;j<m;j++)
		{
			comb[i][j] = 0;
		}
	}
	
	printCombination(comb, arr, ns, m);

	float btemp[m];
	vector < vector<float> > basic_sols(bps, vector<float> (ns, 0)), basic_feasible_sols(bps, vector<float> (n, 0)), temp(m, vector<float> (m, 0));
	for (i=0;i<bps;i++)
	{
		int newtemp[m];
		for (j=0;j<m;j++)
		{
			x = comb[i][j];	
			for (k=0;k<m;k++)
			{							
				temp[k][j] = NewMatrix[k][x];								
			}

			newtemp[j] = x;			
			btemp[j] = b[j];
						
		}
			
		gaussian(m,m,temp,btemp,newtemp,ns,basic_sols,i);
		
	}
	
	int bfps = 0;
	
	for (i=0;i<bps;i++)
	{
		int cfs = 0;
		for (j=0;j<ns;j++)
		{
			if (basic_sols[i][j] > 0)
			{
				cfs++;
			}
		}
		if (cfs == m)
		{
			for (j=0;j<n;j++)
			{	
				basic_feasible_sols[bfps][j] = basic_sols[i][j];				
			}
			bfps++;
		}
	}
	
	float sols[bfps],temps = 0;
	int opt_sols = 0;
	for (i=0;i<bfps;i++)
	{
		for (j=0;j<n;j++)
		{
			temps += basic_feasible_sols[i][j] * Optim[j];
		}
		sols[i] = temps;
		temps = 0;
	}
	
	if (max == 0)
	{
		for (i=0;i<bfps;i++)
		{
			for (j=i;j<bfps;j++)
			{
				if (sols[opt_sols] > sols[j])
				{
					opt_sols = j;
				}
			}
		}
	}
	else if (max == 1)
	{
		for (i=0;i<bfps;i++)
		{
			for (j=i;j<bfps;j++)
			{
				if (sols[opt_sols] < sols[j])
				{
					opt_sols = j;
				}
			}
		}
	}
	else 
	{
		cout << "Problem not properly defined. " << endl;
	}
	
	vector < vector<float> > SimplexTable(m, vector<float> (n+1, 0));
	float Z[n+1], CN[n], CB[m], ratio[m];
	
	string bv[slack_var+1], nbv[n];
	
	CreatSimplexTable(SimplexTable, Z, CN, CB, bv, nbv, n, m, slack_var, NewMatrix, b, Optim, d, a, M);
	
	for (i=0;i<n;i++)
	{
		if (Z[i] < 0)
		{
			Simplex(n, m, SimplexTable, bv, nbv, Z, CN, CB, condn, ratio);
			i=-1;
		}
		if (condn == true)
		{
			break;
		}
	}
	
	while (g != 7){
	switch (g)
	{
		case 1:
			cout << "Basic feasible points are: " << endl;
			for (i=0;i<bfps;i++)
			{
				cout << "( ";
				for (j=0;j<n;j++)
				{
					cout << basic_feasible_sols[i][j] << "  ";
				}
				cout << ")" << endl;
			}
			cout << "\n";
			cout << "Enter next option" << endl;
			cin >> g;
			break;
			
		case 2:
			cout << number_Iteration << endl;
			if (condn == true)
			{
				cout << "Proble of Degeneracy: No optimal solution" << endl;
			}
			cout << "Enter next option" << endl;
			cin >> g;
			break;
			
		case 3:
			int ith;
			cout << "Enter the i value" << endl;
			cin >> ith;
			CreatSimplexTable(SimplexTable, Z, CN, CB, bv, nbv, n, m, slack_var, NewMatrix, b, Optim, d, a, M);
			cout << "Initial Simplex Table: " << endl;
			print_simplex(n, m, SimplexTable, bv, nbv, Z, CN, CB);
			cout << endl;
			j = 0;
			for (i=0;i<n;i++)
			{
				if (j == ith) break;
				
				if (Z[i] < 0)
				{	
					Simplex(n, m, SimplexTable, bv, nbv, Z, CN, CB, condn, ratio);
					j++;
					i=-1;
				}
				if (condn == true)
				{
					break;
				}
			}
			
			cout << "Simplex Table after " << ith << " Iteration: " << endl;
			print_simplex(n, m, SimplexTable, bv, nbv, Z, CN, CB);
			cout << endl;
			
			cout << "Non-basic variables are: " << endl;
			for (i=0;i<n;i++)
			{
				cout << nbv[i];
				cout << " =  " << 0;
				cout << endl;
			}
			cout << "Net evaluation: " << Z[n] << endl;
			cout << "Enter next option" << endl;
			cin >> g;
			break;
    		
		case 4:			
			cout << "Enter the i value" << endl;
			cin >> ith;
			CreatSimplexTable(SimplexTable, Z, CN, CB, bv, nbv, n, m, slack_var, NewMatrix, b, Optim, d, a, M);
			cout << "Initial Simplex Table: " << endl;
			print_simplex(n, m, SimplexTable, bv, nbv, Z, CN, CB);
			cout << endl;
			j = 0;
			for (i=0;i<n;i++)
			{
				if (j == ith) break;
				
				if (Z[i] < 0)
				{	
					Simplex(n, m, SimplexTable, bv, nbv, Z, CN, CB, condn, ratio);
					j++;
					i=-1;
				}
				if (condn == true)
				{
					break;
				}
			}
			
			cout << "Simplex Table after " << ith << " Iteration: " << endl;
			print_simplex(n, m, SimplexTable, bv, nbv, Z, CN, CB);
			cout << endl;
			
			cout << "Basic variables are: " << endl;
			for (i=0;i<slack_var;i++)
			{
				cout << bv[i] << " = " << SimplexTable[i][n] << endl;	
			}
			cout << endl;
			cout << "Ratios are: " << endl;
			for (i=0;i<m;i++)
			{
				cout << ratio[i]  << "   ";
			}
			cout << endl;
			cout << "Enter next option" << endl;
			cin >> g;
			break;
			
		case 5:
			cout << "Enter the i value" << endl;
			cin >> ith;
			cout << "Simplex Table after " << ith << " Iteration: " << endl;
			CreatSimplexTable(SimplexTable, Z, CN, CB, bv, nbv, n, m, slack_var, NewMatrix, b, Optim, d, a, M);
			
			j = 0;
			for (i=0;i<n;i++)
			{
				if (j == ith) break;
				
				if (Z[i] < 0)
				{
					Simplex(n, m, SimplexTable, bv, nbv, Z, CN, CB, condn, ratio);
					j++;
					i=-1;
				}
				if (condn == true)
				{
					break;
				}
			}
			
			print_simplex(n, m, SimplexTable, bv, nbv, Z, CN, CB);
			cout << endl;
			cout << "Enter next option" << endl;
			cin >> g;
			break;
			
		case 6:
			j = 0;
			CreatSimplexTable(SimplexTable, Z, CN, CB, bv, nbv, n, m, slack_var, NewMatrix, b, Optim, d, a, M);
			condn = false;
			for (i=0;i<n;i++)
			{
				if (Z[i] < 0)
				{	
					print_simplex(n, m, SimplexTable, bv, nbv, Z, CN, CB);
					cout << endl;
					Simplex(n, m, SimplexTable, bv, nbv, Z, CN, CB, condn, ratio);
					i=-1;
				}
				if (condn == true)
				{
					cout << endl;
					cout << "Proble of Degeneracy: No optimal solution" << endl;
					cout << endl;
					break;
				}
			}			
			
			if (condn == false)
			{
				cout << "Final Simplex Table is: \n";
				print_simplex(n, m, SimplexTable, bv, nbv, Z, CN, CB);
				cout << endl;
				cout << "Optimal solution: \n[";
				for (i=0;i<m-1;i++)
				{
					cout << bv[i] << ": " << SimplexTable[i][n] << "\n ";
				}
				cout << bv[m-1] << ": " << SimplexTable[m-1][n] << "]" << endl;
			}
			
			for (i=0;i<n;i++)
			{				
				if (Z[i] == 0)
				{	
					alternate = true;
					Simplex(n, m, SimplexTable, bv, nbv, Z, CN, CB, condn, ratio);
				}
				if (condn == true)
				{
					cout << endl;
					cout << "Proble of Degeneracy: No optimal solution" << endl;
					cout << endl;
					break;
				}
			}
			if (alternate == true & condn == false)
			{
				cout << "Alternative Final Simplex Table is: \n";
				print_simplex(n, m, SimplexTable, bv, nbv, Z, CN, CB);
				cout << endl;
				cout << "Alternative Optimal solution: \n[";
				for (i=0;i<m-1;i++)
				{
					cout << bv[i] << ": " << SimplexTable[i][n] << "\n ";
				}
				cout << bv[m-1] << ": " << SimplexTable[m-1][n] << "]" << endl;
			}
			
			cout << "Enter next option" << endl;
			cin >> g;
			break;
		case 7:
			break;
	}}	
	
  		
  	return 0;
}
