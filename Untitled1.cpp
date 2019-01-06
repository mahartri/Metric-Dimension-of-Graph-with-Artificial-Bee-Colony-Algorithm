#include <iostream>
#include <fstream>
#include <string>
using namespace std;

int main()
{
    string Data;
    cout<<" Masukan Nama File "; cin>>Data;
    fstream MyFile;
    MyFile.open(Data, ios::out);
    if(MyFile.is_open())
    {
        cout<<" Performasi Data di Peroleh ";
        cout<<'\n';
    }
    int orde;
    cout<<" Masukan Banyak titik "; cin>>orde;
    cout<<endl;
    int k = 0;
    for(int i = 0; i < orde; i++)
    {
        MyFile<<i+1<<" "<<orde+1<<endl;
    }
    MyFile<<-1<<" "<<-1;
    MyFile.close();
    if(!MyFile.is_open())
    {
        cout<<" File Close ";
        cout<<'\n';
    }


    cout<<" Hello World ";
    system("CLS");
    main();
    return 0;

}
