#include <iostream>
#include"symbol.h"
#include<ctime>

int main() {
    clock_t  begin,end;
    begin=clock();
    computation();
    end=clock();
    std::cout<<"compute over!"<<std::endl;
    std::cout<<"Executed time is "<<(double)(end-begin)/CLOCKS_PER_SEC*1000<<"ms"<<std::endl;
    return 0;
}
