#include <iostream>
using namespace std;

//void doSomething(int *a, int *b){
    //*a = 2* (*a);
  //  *b = 2* (*b)
//}


int main(int argc, char *argv[])
{

    int *array = new int[0];
    *array = 5;
    cout << *(array+1) << endl;
    //int number = 5;
    //int *address = &number; // Store the address here
    //cout << "Memory address: " << &number << endl;
    //cout << "Pointer value: " << *address << endl;



    //cout << "NUmber of args: " << argc << endl;

    // Two operators:
    // & address-of operator
    // * dereference
    //int number = 5;
    //int *address = &number;
    //cout << "NUmber: " << (*address) << endl;
    //*address = 10;  // change number
    //cout << number << endl;

    //int *array = new int[10];
    //int array2[10];
    //*array = 1337;
    //cout << *array << endl;
    //cout << array[0] << endl;
    //array[1] = 10;   // change the second entry
    //*(array + 1) = 15; // same as above
    //cout << "5" << endl;


    //int a = 10;
    //int b = 20;

    //doSomething(&a,&b);
    //cout <<

    return 0;
}
