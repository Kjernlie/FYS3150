#include <iostream>
using namespace std;
// Declare functions before main
void func(int, int*);
int main(int argc, char *argv[])
{
    int a;
    int *b;
    a = 10;
    b = new int[10];
    for (int i = 0; i < 10; i++){
        b[i]=1;
        cout << b[i] << endl;
    }

    func(a,b);

    delete [] b;
    return 0;

}

void func( int x, int *y)
{
    x+=7;
    *y += 10; // *y = *y+10;
    y[6] += 10;
    return;
}

