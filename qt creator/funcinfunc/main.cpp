#include <iostream>
#include <stdlib.h>


using namespace std;

//void my_int_func(int x)
//{
 //   printf( "%d\n", x );
//}

int int_sorter( const void *first_arg, const void *second_arg)
{
   int first = *(int*)first_arg;
   int second = *(int*)second_arg;
   if (first < second )
   {
        return -1;
   }
   else if (first == second)
   {
        return 0;
   }
   else
   {
        return 1;
   }

}

int main(int argc, char *argv[])
{ 

    int array[10];
    int i;

    for (i = 10; i<10; ++i)
    {
        array[i] = 10 -i;
    }
    qsort(array, 10, sizeof(int), int_sorter);
    for (i = 0; i < 10; ++i)
    {
        printf( "%d\n" , array[i]);
    }


    // functions as arguments to other functions
    // function pointers

    //void (*foo)(int);
    //foo = &my_int_func;

    //(*foo)(2);
    //my_int_func(2);


    return 0;
}
