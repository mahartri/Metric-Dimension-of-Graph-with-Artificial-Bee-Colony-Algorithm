#include <windows.h>
#include <iostream>
#include <vector>
#include <cstdlib>
#include <stack>
#include <cmath>
#include <fstream>
#include <ctime>
#include <string>
#include <algorithm>

using namespace std;

//untuk membentuk lebah
struct Bee
{
	vector<double> solusi;
	vector<double> konversi;
	double fungsi_tujuan_sementara;
	double fungsi_tujuan_fix;
	double copy_fungsi_tujuan;
	double trial_limit;
	double fitness;
	double rfitness;
};

//class untuk algoritma ABC
class ArtificialBeeColonyandGraph
{
	public:
		double JumlahIterasi;
		double JumlahSolusi;
		double orde_titik;
		double infinity;
		double jumlah_trial_limit;
		vector<vector<double>> matrix;
		vector<double> Basis;
		vector<Bee> NectarSource;
		vector<Bee> EmployedBee;
		vector<Bee> EmployedBeeBaru;
		vector<Bee> onlookerBee;
		vector<Bee> ScoutBee;
		Bee SaveSolusi;

		ArtificialBeeColonyandGraph(double orde_titik, double JumlahSolusi, double JumlahIterasi, double jumlah_trial_limit);
		//~ArtificialBeeColonyandGraph();
		void Inisialisasi();
		double konversi(double nilai);
		double FungsiTujuanSementera(Bee lebah);
		double FungsiTujuanFix(Bee lebah);
		void CopyFungsiTujuan();
		void EmloyedBeeMencariSolusiNeighbourhood();
		double MencariFitness(Bee lebah);
		void MencariProbabilitas();
		void FaseRoulleteWheel();
		void FaseScoutBee();
		void FaseSaveSolusi();
		void BASIS();
		void addEdge(int start, int end);
		void initializeshortestpath();
		void shortestpath();
		void TampilMatrix();
		void DepthFirstSearch(int Start);
		double random(double start, double end);
};

//inisialisasi solusi awal
ArtificialBeeColonyandGraph::ArtificialBeeColonyandGraph(double orde_titik, double JumlahSolusi, double JumlahIterasi, double jumlah_trial_limit)
{
	this->JumlahSolusi = JumlahSolusi;
	this->orde_titik = orde_titik;
	this->JumlahIterasi = JumlahIterasi;
	this->jumlah_trial_limit = jumlah_trial_limit;

	NectarSource.resize(JumlahSolusi);
	EmployedBee.resize(JumlahSolusi);
	EmployedBeeBaru.resize(JumlahSolusi);
	onlookerBee.resize(JumlahSolusi);
	ScoutBee.resize(JumlahSolusi);
	for(int i = 0; i < JumlahSolusi; i++)
	{
		NectarSource[i].konversi.resize(orde_titik);
		EmployedBee[i].konversi.resize(orde_titik);
		EmployedBeeBaru[i].konversi.resize(orde_titik);
		onlookerBee[i].konversi.resize(orde_titik);
		ScoutBee[i].konversi.resize(orde_titik);
	}
	SaveSolusi.konversi.resize(orde_titik);
	Basis.resize(orde_titik);
	infinity = orde_titik*2;
	matrix.resize(orde_titik,vector<double>(orde_titik));
	for(int i = 0; i < orde_titik; i++)
		{
			for(int j = 0; j < orde_titik; j++)
			{
				matrix[i][j] = 0;
			}
		}
}

//memasukan sisi pada graf
void ArtificialBeeColonyandGraph::addEdge(int start, int end)
{
			if (start > orde_titik || end > orde_titik || start < 0 || end < 0)
			{
				cout<<" invalid edge \n";
			}
			else if(start == end)
			{
				matrix[start-1][end-1] = matrix[start-1][end-1]= 0;
			}
			else
			{
				matrix[start-1][end-1] = matrix[end-1][start-1]= 1;
			}
}

//Algoritma mengujungi titik titik pada graf
void ArtificialBeeColonyandGraph::DepthFirstSearch (int start)
{
		stack<int> s;
		vector<bool> visited;
		vector<bool> visit;

		for (int i = 0; i < orde_titik; i++)
		{
			visited.push_back(false);
			visit.push_back(false);
		}

		s.push(start);
		visited[start-1] = true;
		cout << " DFS : " << start<<" ";
		while (!s.empty())
		{
			for (int i = 0; i < orde_titik; i++)
			{
				if ((matrix[start-1][i] != 0) && (visited[i] == false) && (visit[i]== false))
				{
					s.push(i+1);
					visit[i] = true;
					start = i+1;
					cout<<start<<" ";

				}
			}
			start = s.top();
			s.pop();
			visit[start-1] = false;
			visited[start-1] = true;
		}
}

//menginisialisasi matriks untuk menentukan lintasan terpendek
void ArtificialBeeColonyandGraph::initializeshortestpath()
{
	for(int i = 0 ; i < orde_titik; i++)
		{
			for(int j = 0; j < orde_titik; j++)
			{
				if((i!=j) && (matrix[i][j] == 0))
				{
					matrix[i][j] = infinity;
				}
			}
		}
}

//mencari lintasan terpendek
void ArtificialBeeColonyandGraph::shortestpath()
{
	for(int k = 0; k < orde_titik; k++)
		{
			for(int i = 0; i < orde_titik; i++)
			{
				for(int j = 0; j < orde_titik; j++)
				{
					if(matrix[i][k]+matrix[k][j] < matrix[i][j])
					{
						matrix[i][j] = matrix[i][k]+matrix[k][j];
					}
				}
			}
		}
}

//untuk menampilkan hasil lintasan terpendek
void ArtificialBeeColonyandGraph::TampilMatrix()
{
		for(int i = 0; i < orde_titik; i++)
		{
			for(int j = 0; j < orde_titik; j++)
			{
				cout<<matrix[i][j]<<"\t";
			}
			cout<<"\n";
		}
}

//untuk melakukan random solusi dengan [start, end]
double ArtificialBeeColonyandGraph::random(double start, double end)
{
	return start + (end-start)*rand()/(RAND_MAX + 1.0);
}

//untuk mengkonversi solusi ke angka 0 dan 1
double ArtificialBeeColonyandGraph::konversi(double nilai)
{
	double konversi1;
		if((nilai <= 0.5) && (nilai >= -0.5))
			{
				konversi1 = 1;
			}

		else if((nilai > 0.5) || (nilai < -0.5))
			{
				konversi1 = 0;
			}
	return konversi1;
}

//menentukan fungsi tujuan sementara
double ArtificialBeeColonyandGraph::FungsiTujuanSementera(Bee lebah)
{
    int sum = 0;
    for(int j = 0; j < orde_titik; j++)
        {
            if(lebah.konversi[j] == 1)
                sum++;
        }
  	return sum;
}


//menentukan fungsi tujuan yang sebenarnya
double ArtificialBeeColonyandGraph::FungsiTujuanFix(Bee lebah)
{
	vector<vector<int>> representasi;
	for (int i = 0; i < orde_titik; i++)
	{
		vector<int> temp;
		for( int j = 0; j < orde_titik; j++)
		{
			if(lebah.konversi[j] == 1)
			{
				temp.push_back(matrix[i][j]);
			}
		}
		representasi.push_back(temp);
	}

	int fungsitujuanfix = 0;
	int uniksemua = 0;
	int total_perbandingan = 0;

	for(int i = 0; i < orde_titik; i++)
	{
		int k = i+1;
		while(k < orde_titik)
		{
			total_perbandingan++;
			int loop = 0;
			int parameter_kembar = 0;
			while( loop < lebah.fungsi_tujuan_sementara)
			{
				if(representasi[i][loop] != representasi[k][loop] )
				{
					parameter_kembar = 0;
					uniksemua++;
					loop = lebah.fungsi_tujuan_sementara+1;
					k++;
				}
				else
				{
					parameter_kembar++;
					loop++;
				}
			}
			if(lebah.fungsi_tujuan_sementara == parameter_kembar)
			{
				fungsitujuanfix = orde_titik + 1;
				k = orde_titik + 1;
				i = orde_titik + 1;
			}
		}
	}
	if( uniksemua == total_perbandingan)
		{
			fungsitujuanfix = lebah.fungsi_tujuan_sementara;
		}
	return fungsitujuanfix;
}

void ArtificialBeeColonyandGraph::BASIS()
{
	vector<double> basis1;
	for(int j = 0; j < orde_titik; j++)
            {
                if(SaveSolusi.konversi[j] == 1)
                {
                    basis1.push_back(j+1);
                }
            }

    if(basis1.size() > orde_titik)
    	basis1.pop_back();

    cout<<endl;
    cout<<basis1.size();
    cout<<" ( ";
    for(int i = 0; i < SaveSolusi.fungsi_tujuan_fix; i++)
    {
    	cout<<basis1[i]<<" ";
    }
    cout<<" ) ";
}

//menginisialisasi parameter awal ABC
void ArtificialBeeColonyandGraph::Inisialisasi()
{
	for(int i = 0; i < JumlahSolusi; i++)
	{
		for(int j = 0; j < orde_titik; j++)
		{
			NectarSource[i].solusi.push_back(random(-1, 1));
			EmployedBee[i].solusi.push_back(NectarSource[i].solusi[j]);
			EmployedBeeBaru[i].solusi.push_back(NectarSource[i].solusi[j]);
			onlookerBee[i].solusi.push_back(NectarSource[i].solusi[j]);
			ScoutBee[i].solusi.push_back(NectarSource[i].solusi[j]);
			SaveSolusi.solusi.push_back(NectarSource[0].solusi[j]);
		}
		//inisialisasi fungsi tujuan
		for(int j = 0; j < orde_titik; j++)
		{
			NectarSource[i].konversi[j] = konversi(NectarSource[i].solusi[j]);
		}
		NectarSource[i].fungsi_tujuan_sementara = FungsiTujuanSementera(NectarSource[i]);
		NectarSource[i].fungsi_tujuan_fix = FungsiTujuanFix(NectarSource[i]);
		NectarSource[i].fitness = MencariFitness(NectarSource[i]);
		NectarSource[i].rfitness = 0;
		NectarSource[i].trial_limit = 0;

		//inisialisasi Employed Bee
		for(int j = 0; j < orde_titik; j++)
		{
			EmployedBee[i].konversi[j] = NectarSource[i].konversi[j];
		}
		EmployedBee[i].fungsi_tujuan_sementara = NectarSource[i].fungsi_tujuan_sementara;
		EmployedBee[i].fungsi_tujuan_fix = NectarSource[i].fungsi_tujuan_fix;
		EmployedBee[i].trial_limit = NectarSource[i].trial_limit;
		EmployedBee[i].fitness = NectarSource[i].fitness;
		EmployedBee[i].rfitness = NectarSource[i].rfitness;

		//Inisialisasi Onlooker Bee
		for(int j = 0; j < orde_titik; j++)
		{
			onlookerBee[i].konversi[j] = NectarSource[i].konversi[j];
		}
		onlookerBee[i].fungsi_tujuan_sementara = NectarSource[i].fungsi_tujuan_sementara;
		onlookerBee[i].fungsi_tujuan_fix = NectarSource[i].fungsi_tujuan_fix;
		onlookerBee[i].trial_limit = NectarSource[i].trial_limit;
		onlookerBee[i].fitness = NectarSource[i].fitness;
		onlookerBee[i].rfitness = NectarSource[i].rfitness;

		//Inisialisasi Scout Bee
		for(int j = 0; j < orde_titik; j++)
		{
			ScoutBee[i].konversi[j] = NectarSource[i].konversi[j];
		}
		ScoutBee[i].fungsi_tujuan_sementara = NectarSource[i].fungsi_tujuan_sementara;
		ScoutBee[i].fungsi_tujuan_fix = NectarSource[i].fungsi_tujuan_fix;
		ScoutBee[i].trial_limit = NectarSource[i].trial_limit;
		ScoutBee[i].fitness = NectarSource[i].fitness;
		ScoutBee[i].rfitness = NectarSource[i].rfitness;
	}
	for(int j = 0; j < orde_titik; j++)
		{
			SaveSolusi.konversi[j] = NectarSource[0].konversi[j];
		}
	SaveSolusi.fungsi_tujuan_sementara = NectarSource[0].fungsi_tujuan_sementara;
	SaveSolusi.fungsi_tujuan_fix = NectarSource[0].fungsi_tujuan_fix;
	SaveSolusi.fitness = NectarSource[0].fitness;
	SaveSolusi.rfitness = NectarSource[0].rfitness;
	SaveSolusi.trial_limit = NectarSource[0].trial_limit;
}

//Mencari Persekitaran dari Solusi Lama
void ArtificialBeeColonyandGraph::EmloyedBeeMencariSolusiNeighbourhood()
{
	int k;
	int titikyangingindiganti;
	double nilairandomij; //nilai random antara [-1, 1]
	for(int i = 0; i < JumlahSolusi; i++)
	{
		for(int j = 0; j < orde_titik; j++)
			{
				EmployedBee[i].solusi[j] = NectarSource[i].solusi[j];
			}

		for(int j = 0; j < orde_titik; j++)
		{
			//memilih suatu solusi lama yang akan digunakan untuk mencari persekitaran dari solusi ke i;
			//solusi lama yang terpilih secara acak haruslah berbeda dengan solusi ke - i
			while(true)
			{
				k = (int)random(0, JumlahSolusi);
				if(k != i)
					break;
			}
			//memilih titik secara acak pada solusi ke - i yang akan diganti nilai nya
			titikyangingindiganti = (int)random(0, orde_titik); //"(int)" mengkonversi nilai dari suatu tipe data ke tipe data "int"
			nilairandomij = random(-1, 1);

			EmployedBee[i].solusi[titikyangingindiganti] = NectarSource[i].solusi[titikyangingindiganti] + nilairandomij*(NectarSource[i].solusi[titikyangingindiganti] - NectarSource[k].solusi[titikyangingindiganti]);

			if (EmployedBee[i].solusi[titikyangingindiganti] < -1)
			{
				EmployedBee[i].solusi[titikyangingindiganti] = -1;
			}
			if (EmployedBee[i].solusi[titikyangingindiganti] > 1)
			{
				EmployedBee[i].solusi[titikyangingindiganti] = 1;
			}
		}

		//menghitung nilai fungsi tujuan sementara dan nilai fungsi tujuan baru
		for(int j = 0; j < orde_titik; j++)
		{
			EmployedBee[i].konversi[j] = konversi(EmployedBee[i].solusi[j]);
		}
		EmployedBee[i].fungsi_tujuan_sementara = FungsiTujuanSementera(EmployedBee[i]);
		EmployedBee[i].fungsi_tujuan_fix = FungsiTujuanFix(EmployedBee[i]);
		EmployedBee[i].fitness = MencariFitness(EmployedBee[i]);

		// uji perbandingan solusi baru (neighbourhood) dengan solusi lama
		// jika solusi baru lebih baik dari solusi lama maka solusi lama diganti dengan solusi baru dan triall limit solusi baru di reset ke nol
		// jika tidak, maka solusi lama dipertahankan dan trial limit solusi lama ditambah 1
		if (EmployedBee[i].fungsi_tujuan_fix < NectarSource[i].fungsi_tujuan_fix)
		{
			for(int j = 0 ; j < orde_titik; j++)
			{
				NectarSource[i].solusi[j] = EmployedBee[i].solusi[j];
				NectarSource[i].konversi[j] = EmployedBee[i].konversi[j];
			}
			NectarSource[i].trial_limit = 0;
			NectarSource[i].fungsi_tujuan_sementara = EmployedBee[i].fungsi_tujuan_sementara;
			NectarSource[i].fungsi_tujuan_fix = EmployedBee[i].fungsi_tujuan_fix;
		}
		else
		{
			NectarSource[i].trial_limit++;
		}
	}
}

//Mencari nilai fitnes dari solusi ke i
double ArtificialBeeColonyandGraph::MencariFitness(Bee lebah)
{
	double temp;
		temp = 1/lebah.fungsi_tujuan_fix;
		return temp;
}

void ArtificialBeeColonyandGraph::MencariProbabilitas()
{
	//mencari nilai total fungsi tujuan
	//double ZigmaNilaiFitnes;
	//ZigmaNilaiFitnes = NectarSource[0].fitness;
	//for(int i = 1; i < JumlahSolusi; i++)
		//ZigmaNilaiFitnes += NectarSource[i].fitness;

	double maxfitness;
	maxfitness = NectarSource[0].fitness;
	for(int i = 1; i < JumlahSolusi; i++)
	{
	if(NectarSource[i].fitness > maxfitness)
			maxfitness = NectarSource[i].fitness;
	}

	for( int i = 0; i < JumlahSolusi; i++)
		NectarSource[i].rfitness = (0.9*(NectarSource[i].fitness/maxfitness)) + 0.1;
	//for(int i = 0 ; i < JumlahSolusi; i++)
		//NectarSource[i].rfitness = NectarSource[i].fitness/ZigmaNilaiFitnes;
}

void ArtificialBeeColonyandGraph::FaseRoulleteWheel()
{
	int k, i, t;
	i = 0;
	t = 0;

	double rfitnessDipilih;
	int titikyangingindiganti;
	double nilairandomij; //nilai random antara [-1, 1]

	while(t < JumlahSolusi)
	{
		rfitnessDipilih = random(0, 1);
		//memilih solusi yang merupakan domain dari rfitnessDipilih
		if(rfitnessDipilih < NectarSource[i].rfitness)
		{
			t++;
			for(int j = 0; j < orde_titik; j++)
				{
					onlookerBee[i].solusi[j] = NectarSource[i].solusi[j];
				}

			for(int j = 0; j < orde_titik; j++)
			{
				//memilih suatu solusi lama yang akan digunakan untuk mencari persekitaran dari solusi ke i;
				//solusi lama yang terpilih secara acak haruslah berbeda dengan solusi ke - i
				while(true)
				{
					k = (int)random(0, JumlahSolusi);
					if(k != i)
						break;
				}
				//memilih titik secara acak pada solusi ke - i yang akan diganti nilai nya
				titikyangingindiganti = (int)random(0, orde_titik); //"(int)" mengkonversi nilai dari suatu tipe data ke tipe data "int"
				nilairandomij = random(-1, 1);

				onlookerBee[i].solusi[titikyangingindiganti] = NectarSource[i].solusi[titikyangingindiganti] + nilairandomij*(NectarSource[i].solusi[titikyangingindiganti] - NectarSource[k].solusi[titikyangingindiganti]);

				if (onlookerBee[i].solusi[titikyangingindiganti] < -1)
				{
					onlookerBee[i].solusi[titikyangingindiganti] = -1;
				}
				if (onlookerBee[i].solusi[titikyangingindiganti] > 1)
				{
					onlookerBee[i].solusi[titikyangingindiganti] = 1;
				}
			}

			//menghitung nilai fungsi tujuan sementara dan nilai fungsi tujuan baru
			for(int j = 0; j < orde_titik; j++)
			{
				onlookerBee[i].konversi[j] = konversi(onlookerBee[i].solusi[j]);
			}
			onlookerBee[i].fungsi_tujuan_sementara = FungsiTujuanSementera(onlookerBee[i]);
			onlookerBee[i].fungsi_tujuan_fix = FungsiTujuanFix(onlookerBee[i]);
			onlookerBee[i].fitness = MencariFitness(onlookerBee[i]);

			// uji perbandingan solusi baru (neighbourhood) dengan solusi lama
			// jika solusi baru lebih baik dari solusi lama maka solusi lama diganti dengan solusi baru dan triall limit solusi baru di reset ke nol
			// jika tidak, maka solusi lama dipertahankan dan trial limit solusi lama ditambah 1
			if (onlookerBee[i].fungsi_tujuan_fix < NectarSource[i].fungsi_tujuan_fix)
			{
				for(int j = 0 ; j < orde_titik; j++)
				{
					NectarSource[i].solusi[j] = onlookerBee[i].solusi[j];
					NectarSource[i].konversi[j] = onlookerBee[i].konversi[j];
				}
				NectarSource[i].trial_limit = 0;
				NectarSource[i].fungsi_tujuan_sementara = onlookerBee[i].fungsi_tujuan_sementara;
				NectarSource[i].fungsi_tujuan_fix = onlookerBee[i].fungsi_tujuan_fix ;
				NectarSource[i].fitness = onlookerBee[i].fitness;
			}
			else
			{
				NectarSource[i].trial_limit++;
			}
		}
		i++;
		if(i == JumlahSolusi)
			i = 0;
	}
}

//scoutbee mencari solusi baru
void ArtificialBeeColonyandGraph::FaseScoutBee()
{
	int Temp1;
	Temp1 = 0;
	for(int i = 0 ; i < JumlahSolusi; i++)
	{
		if(NectarSource[i].trial_limit > NectarSource[Temp1].trial_limit)
			Temp1 = i;
	}

	//membangun soulsi baru pada solusi lama yang memiliki trial limit >= limit
	if(NectarSource[Temp1].trial_limit >= jumlah_trial_limit)
	{
		for(int i = 0; i < orde_titik; i++)
		{
			SaveSolusi.solusi[i] = NectarSource[Temp1].solusi[i];
			SaveSolusi.konversi[i] = NectarSource[Temp1].konversi[i];
		}
		SaveSolusi.fungsi_tujuan_sementara =  NectarSource[Temp1].fungsi_tujuan_sementara;
		SaveSolusi.fungsi_tujuan_fix = NectarSource[Temp1].fungsi_tujuan_fix;
		SaveSolusi.fitness = NectarSource[Temp1].fitness;

		for(int i = 0; i < orde_titik; i++)
		{
			double RANDOM;
			RANDOM = random(-1,1);
			NectarSource[Temp1].solusi[i] = RANDOM;
		}
		NectarSource[Temp1].trial_limit = 0;
		for(int j = 0; j < orde_titik; j++)
		{
			NectarSource[Temp1].konversi[j] = konversi(NectarSource[Temp1].solusi[j]);
		}
		NectarSource[Temp1].fungsi_tujuan_sementara = FungsiTujuanSementera(NectarSource[Temp1]);
		NectarSource[Temp1].fungsi_tujuan_fix = FungsiTujuanFix(NectarSource[Temp1]);
		NectarSource[Temp1].fitness = MencariFitness(NectarSource[Temp1]);
	}
	else
	{
		MencariProbabilitas();
		FaseRoulleteWheel();
		FaseScoutBee();
	}

}

//menyimpan solusi terbaik sementara
void ArtificialBeeColonyandGraph::FaseSaveSolusi()
{
	for (int i = 0; i < JumlahSolusi; i++)
	{
		if (NectarSource[i].fungsi_tujuan_fix < SaveSolusi.fungsi_tujuan_fix)
		{
			for (int j = 0; j<orde_titik; j++)
			{
				SaveSolusi.solusi[j] = NectarSource[i].solusi[j];
				SaveSolusi.konversi[j] = NectarSource[i].konversi[j];
			}
			SaveSolusi.fungsi_tujuan_fix = NectarSource[i].fungsi_tujuan_fix;
		}
	}
}

int main(void)
{
	clock_t tStart = clock();
	cout<<""<<endl;
	cout<<"PENENTUAN DIMENSI METRIK GRAF MENGGUNAKAN ALGORITMA ARTIFICIAL BEE COLONY";
	int awal;
    double BanyakTrialLimit, BanyakIterasi;
    cout << endl;

    double BanyakSolusi;
    cout<<"\nMasukan jumlah solusi 				= ";
    cin>>BanyakSolusi;
    if(BanyakSolusi < 1)
    {
    	cout<<" \n input yang anda masukan salah \n";
    	system("PAUSE");
    	system("CLS");
    	main();
    }

    cout<<"\nMasukan jumlah banyak iterasi 			= "; cin>>BanyakIterasi;
    if(BanyakIterasi < 1)
    {
    	cout<<" \n input yang anda masukan salah \n";
    	system("PAUSE");
    	system("CLS");
    	main();
    }

    cout<<"\nmasukan jumlah banyak trial imit 		= "; cin>>BanyakTrialLimit;
    if(BanyakTrialLimit < 1)
    {
    	cout<<" \n input yang anda masukan salah \n";
    	system("PAUSE");
    	system("CLS");
    	main();
    }

    int asal,tujuan;
   	string dirGraf;
    cout<<" Masukan Directory(/alamat/nama) file Graf : ";
    cin>>dirGraf;
    cout<<endl;

    vector<double> tampung;
    ifstream file (dirGraf);
    if(file.is_open())
    	cout<<" File has opened \n";

    while(!file.eof())
    {
    	file>>asal>>tujuan;
    	if((asal == -1) && (tujuan == -1))
        {
            break;
        }
    	tampung.push_back(asal);
        tampung.push_back(tujuan);
    }

    double MaxElement1 = tampung[0];
    for(int i = 1; i < tampung.size(); i++)
    {
    	if(MaxElement1 < tampung[i])
    		MaxElement1 = tampung[i];
    }
    cout<<"\nMasukan jumlah order titik 			= ";
    cout<<MaxElement1;

    ArtificialBeeColonyandGraph LEBAH(MaxElement1,BanyakSolusi, BanyakIterasi, BanyakTrialLimit);
	cout<<endl;
    ifstream file1 (dirGraf);
    while(!file1.eof())
    {
        cout<<" enter edge (-1, -1 to exit) : ";
        file1>>asal>>tujuan;
        if((asal == -1) && (tujuan == -1))
        {
            break;
        }
        cout<<asal<<" "<<tujuan;
        LEBAH.addEdge(asal, tujuan);
        cout<<endl;
    }


    cout<<"\n\n";
    LEBAH.TampilMatrix();
    cout<<endl;
    cout<<" masukan titik awal yang dikunjungi ";
    cin>>awal;
    cout<<"\n\n";
    LEBAH.DepthFirstSearch(awal);
    LEBAH.initializeshortestpath();
    LEBAH.shortestpath();
    LEBAH.TampilMatrix();
    double noiterasi;
    noiterasi = 0;
    LEBAH.Inisialisasi();
    for(int i = 0 ; i < BanyakSolusi; i++)
    {
    	for(int j = 0; j < MaxElement1; j++)
    	{
    		cout<<LEBAH.EmployedBee[i].konversi[j]<<" ";
    	}
    	cout<<endl;
    }
    cout<<endl;
  	for(int i = 0; i < BanyakSolusi; i++)
    {
    	for(int j = 0; j < MaxElement1; j++)
    	{
    		cout<<LEBAH.EmployedBee[i].solusi[j]<<" ";
    	}
    	cout<<endl;
    	for(int j = 0; j < MaxElement1; j++)
    	{
    		cout<<LEBAH.EmployedBee[i].konversi[j]<<" ";
    	}
    	cout<<endl;
    	for(int j = 0; j < MaxElement1; j++)
    	{
    		cout<<LEBAH.NectarSource[i].konversi[j]<<" ";
    	}
    	cout<<endl;
    	cout<<LEBAH.EmployedBee[i].fungsi_tujuan_fix<<" ";
    }
    cout<<LEBAH.NectarSource[0].fungsi_tujuan_sementara;
    LEBAH.FaseSaveSolusi();
    vector<vector<double>> SaveKonversi;
    vector<float> SaveFungsiTujuan;
    SaveFungsiTujuan.resize(BanyakIterasi);
    while(noiterasi < BanyakIterasi)
    {
    	vector<double> temp;
    	LEBAH.EmloyedBeeMencariSolusiNeighbourhood();

		LEBAH.MencariProbabilitas();

		LEBAH.FaseRoulleteWheel();

		LEBAH.FaseSaveSolusi();

		LEBAH.FaseScoutBee();

		SaveFungsiTujuan[noiterasi] = LEBAH.SaveSolusi.fungsi_tujuan_fix;
		for(int i = 0; i < MaxElement1; i++)
		{
			temp.push_back(LEBAH.SaveSolusi.konversi[i]);
		}

		for(int i = 0; i < MaxElement1; i++)
		{
			cout<<temp[i]<<" ";
		}
		SaveKonversi.push_back(temp);
		cout<<" \nItersai ke-"<<noiterasi;
		cout<<" \nDimensi Metrik sementara yang diperoleh 				= "<<LEBAH.SaveSolusi.fungsi_tujuan_fix;
		cout<<endl;
		noiterasi++;
    }
    cout<<endl;

    int temp = SaveFungsiTujuan[0];
    vector<double> tempkonversi;
    tempkonversi.resize(MaxElement1);
    for(int i = 1; i < BanyakIterasi; i++)
    {
    	vector<double> temporary;
    	if( temp > SaveFungsiTujuan[i])
    		{
    			temp = SaveFungsiTujuan[i];
    			for(int j = 0; j < MaxElement1; j++)
    			{
    				SaveKonversi[0][j] = SaveKonversi[i][j];
    			}
    		}
    }

    cout<<" Diperoleh dimensi metrik-nya adalah : "<<temp;
    cout<<" \nDan basisnya adalah "<<" ( ";
    for(int i = 0; i < MaxElement1; i++)
    {
    	if(SaveKonversi[0][i] == 1)
    		cout<<i+1<<" ";
    }
 	cout<<" ) ";
    cout<<endl;
	cout<<"hello world ! "<<endl;
	cout<<"Time taken: \n"<<(double)(clock() - tStart)/CLOCKS_PER_SEC;
	cout<<endl;
	system("PAUSE");
	cin.get();
	system("CLS");
	main();
	return 0;
}
