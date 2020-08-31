#include <iostream>
#include <vector>
#include<cstring>
#include <string>
#include <malloc.h>
//#include <malloc_np.h>
#include<stdlib.h>
using namespace std;

template<class T>
int _msize(T *p){

}

int main()
{
    vector<string> msg {"Hello", "C++", "World", "from", "VS Code", "and the C++ extension!"};

	int qq = 21;
	double **p = new double *[qq];
	for (int i = 0; i < qq; i++)
	{
		p[i] = new double[40];
		memset(p[i], 0, 5 * sizeof(double));
	}
	cout << malloc_usable_size(p)/sizeof(p[0]) << endl;
    cout << (p[malloc_usable_size(p)/sizeof(p)-1] == p[malloc_usable_size(p)/sizeof(p)-2]) << endl;
	cout << malloc_usable_size(p[0])/sizeof(p[0][0]) << endl;
	cout << p[21] << endl;
}