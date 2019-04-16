#include <iostream>
#include <fstream>
#include <string>
#include <vector>

using namespace std;
#define LONG_LONG_BITS 64

// 64-bit random seeds corresponding to bases and their complements
static const uint64_t seedA = 0x3c8bfbb395c60474;
static const uint64_t seedC = 0x3193c18562a02b4c;
static const uint64_t seedG = 0x20323ed082572324;
static const uint64_t seedT = 0x295549f54be24456;

uint64_t primo_Hash(string sequenza_input, int spaceSeedSize);
uint64_t restanti_Hash(uint64_t precedente, char uscente, char entrante, int spaceSeedSize);
vector <uint64_t> total_seq_hash(string sequenza_input, int spaceSeedSize);
string complemento(string spacedSeed);
vector<int> preProcessing_1(string spacedSeedComp);
vector<int> preProcessing_2(string spacedSeedComp, vector<int> ab);
uint64_t toInt(char seq);
uint64_t leftRotate(uint64_t n, int d);
uint64_t rightRotate(uint64_t n, int d);
vector<string> read_save_file(string nomefile);

int main()
{
	//parse command line DA FARE//////////
	string nomefile = "C:\\Users\\ricca\\Desktop\\2° anno 2° semestre\\BioInfo\\single_end_dataset\\R1.fna";
	//apre file 
	ifstream inputfile;
	inputfile.open(nomefile);
	vector<string> file = read_save_file(nomefile);


	string spacedSeed = "1111111111111111111010101";
	int spaceSeedSize = spacedSeed.length();

	for (int i = 0; i < file.size(); i++) {
		cout << file[i] << endl;
	}

	/*
	vector<uint64_t> hashVector;
	string prova = "ACGTACGTACTGTCAGTCGATGCCTAGGTCATGCATGCGTCATGCATCATGCAGTCAGTACGTGCTATCGATGACTGACCAGTCAGTCGTAGCTAGCTAGCTAGCTCATGACTACTGACTATCTACTACTACTATCACC";
	hashVector = total_seq_hash(prova, spaceSeedSize);
	for (int i = 0; i < hashVector.size(); i++) {
		cout << hashVector[i] << endl;
	}
	*/
	/*
	string comp = complemento(spacedSeed);
	cout << comp << endl;
	cout << "preprocessing 2  " << endl;
	vector<int> p;
	p = preProcessing_1(comp);
	vector<int> s;
	s = preProcessing_2(comp,p);
	for (int i = 0; i < s.size(); i++)
	{
		cout << s[i];
		cout << " ";
	}

	*/


	cin.ignore();
}

//trova il primo hash, dalla posizione 0 a k-1
uint64_t primo_Hash(string sequenza_input, int spaceSeedSize) {

	//il primo carattere
	uint64_t hash = leftRotate(toInt(sequenza_input[0]), spaceSeedSize - 1);
	uint64_t temp;
	//ciclo lungo come il k-mer
	for (int i = 1; i < spaceSeedSize; i++)
	{
		//calcola lo shift necessario e poi fa lo XOR con quanto già calcolato 
		temp = leftRotate(toInt(sequenza_input[i]), spaceSeedSize - i + 1);
		hash = hash ^ temp;
	}

	return hash;
}

//trova i restanti hash fino all'ultimo k-mer
uint64_t restanti_Hash(uint64_t precedente, char uscente, char entrante, int spaceSeedSize) {

	//hash di quello che sta prima ruotato di 1 XOR hash carattere che esce ruotato di k XOR hash di quello che esce 
	uint64_t hash = leftRotate(precedente, 1) ^ leftRotate(toInt(uscente), spaceSeedSize) ^ toInt(entrante);
	return hash;
}

//trova tutti gli hash e restituisce il vettore con gli hash, la posizione è l'indice di hash
vector <uint64_t> total_seq_hash(string sequenza_input, int spaceSeedSize) {

	//dentro il vettore ci salvo tutti gli hash 
	vector <uint64_t> hashVector(sequenza_input.length() - spaceSeedSize + 1);
	//per fare il primo hash mi servono solo k valori, taglio la sequenza di conseguenza
	string primicaratteri = sequenza_input.substr(0, spaceSeedSize - 1);
	//calcolo l'hash del primo elemento e lo salvo
	uint64_t primoHash = primo_Hash(primicaratteri, spaceSeedSize);
	hashVector[0] = primoHash;

	//calcolo tutti gli altri hash e li aggiungo nel vettore
	//indice parte da 1 perchè 0 già calcolato,
	//si ferma ad i = (l-k) (lunghezza input- lunghezza finestra)
	for (int i = 1; i <= sequenza_input.length() - spaceSeedSize; i++)
	{
		//hash di quello che sta prima, il primo carattere che esce, il primo carattere che entra nella finestra lunga k  
		hashVector[i] = restanti_Hash(hashVector[i - 1], sequenza_input[i], sequenza_input[spaceSeedSize + i - 1], spaceSeedSize);
	}

	return hashVector;
}

//dalla sequenza prendere il relativo "hash numerico"
uint64_t toInt(char seq) {
	if (seq == 'A') return seedA;
	if (seq == 'C') return seedC;
	if (seq == 'G') return seedG;
	if (seq == 'T') return seedT;
	//errore
	return -1;
}

//fa la rotazione di n per d bits
uint64_t leftRotate(uint64_t n, int d)
{
	return (n << d) | (n >> (LONG_LONG_BITS - d));
}

//idem solo a destra
uint64_t rightRotate(uint64_t n, int d)
{
	return (n >> d) | (n << (LONG_LONG_BITS - d));
}

/*
//calcola la stringa complemento della stringa spaceSeed
uint64_t hash_complemento(string spacedSeed, string sequenza_input)
{
	//calcolo il complemento che mi serve
	int pos1 = 0;
	int pos2 = 0;
	bool primo = true;
	for	(int i = 0; i < spacedSeed.length(); i++) {
		//cerco la pos del primo zero
		if (primo) {
			if (spacedSeed[i] == '0') {
				pos1 = i;
				primo = false;
			}

		}
		if (spacedSeed[i] == '0')
			pos2 = i;
	}
	//error
	if (pos1 == 0 && pos2 == 0) return -1;

	//creo la stringa
	string complemento;
	for (int i = pos1; i <= pos2; i++) {
		if (spacedSeed[i] == '0')
			complemento.append("1");
		if (spacedSeed[i] == '1')
			complemento.append("0");
	}

	//calcolo l'hash della stringa complemento
	uint64_t hash_complemento = 4445;
	return hash_complemento;
}
*/
/*
////finding the reduced complementary of the spaced seed
string complementary(string spacedSeed) {
	string comp;
	for (int i = 0; i < spacedSeed.length(); i++) {
		if (spacedSeed[i] == '0')
			comp.push_back('1');
		else
			comp.push_back('0');
	}

	int i = 0;
	int j = spacedSeed.length() - 1;
	while (comp[i] == '0') {
		i++;
	}
	while (comp[j] == '0') {
		j--;
	}
	comp.erase(0, i);
	comp.erase(j - (i - 1), comp.length());
	return comp;
}
*/

//ritorn la stinga complementata con le x all'inizio e alla fine
string complemento(string spacedSeed)
{
	//indice degli zeri 
	int primo;
	int ultimo;
	bool x = true;

	//cerrco primo e ultimo
	for (int i = 0; i < spacedSeed.length(); i++) {
		if (spacedSeed[i] == '0' && x) {
			primo = i;
			x = false;
		}
		if (spacedSeed[i] == '0' && (!x)) {
			ultimo = i;
		}
	}

	//metto le x all'inizio e alla fine 
	for (int i = 0; i < primo; i++)
		spacedSeed[i] = 'x';
	for (int i = ultimo + 1; i < spacedSeed.length(); i++)
		spacedSeed[i] = 'x';

	//agli altri faccio girare 
	for (int i = primo; i <= ultimo; i++) {
		if (spacedSeed[i] == '0')
			spacedSeed[i] = '1';
		else
			spacedSeed[i] = '0';
	}

	return spacedSeed;
}

////preprocessing numero 1 per trovare a e b
vector<int> preProcessing_1(string spacedSeedComp) {

	//posizione dove conviene attaccarsi
	vector<int> ab(2);
	int score = 0;
	int bestScore = score;

	int primo;
	int ultimo;
	bool x = true;

	//cerrco primo e ultimo  x 
	for (int i = 0; i < spacedSeedComp.length(); i++) {
		if (spacedSeedComp[i] == 'x' && x) {
			primo = i;
		}
		else {
			x = false;
		}
		if (spacedSeedComp[i] == 'x' && x != true) {
			ultimo = i;
			break;
		}
	}

	ab[1] = primo + 1;
	//creo string per fare i confronti
	string temp;
	for (int i = primo + 1; i < ultimo; i++)
		temp.push_back(spacedSeedComp[i]);

	//numero di sovrapposizioni possibili (ultimo -primo+2)
	for (int i = primo + 2; i < ultimo; i++) {

		for (int j = 0; j < temp.length(); j++) {
			//passato la soglia 
			if ((i + j) >= ultimo) break;
			//conmti somiglianze
			if (spacedSeedComp[j + i] == temp[j])
				score++;
		}
		if (score > bestScore) {
			bestScore = score;
			//punto di attacco del primo 
			ab[0] = i;
		}
		score = 0;
	}

	return ab;

}

//secondo preprocessing per fare il vettore di comandi, 0 niente, 1 ins, 2 canc, 3 ins/canc
vector<int> preProcessing_2(string spacedSeedComp, vector<int> ab) {

	vector<int> comandi(spacedSeedComp.length());
	//creo le due stringhe in modo da capire che comandi assegnare
	string primo = spacedSeedComp.substr(ab[0], spacedSeedComp.length() - ab[0]) + spacedSeedComp.substr(0, ab[0]);
	string secondo = spacedSeedComp.substr(ab[1], spacedSeedComp.length() - ab[1]) + spacedSeedComp.substr(0, ab[1]);

	//la prima parte se è gia apposto non si tocca il resto si
	for (int i = 0; i < primo.length(); i++) {
		if (primo[i] == 'x')
			break;
		if ((primo[i] == '1') && (secondo[i] == '1')) {
			primo[i] = '0';
			secondo[i] = '0';
		}
	}

	//la seconda parte si aggiusta di conseguenza mettendo il comando corretto
	for (int i = 0; i < primo.length(); i++) {
		//caso (0,0) (0,x) (x,0) (x,x) 
		if (((primo[i] == 'x') || (primo[i] == '0')) && ((secondo[i] == 'x') || (secondo[i] == '0'))) {
			comandi[i] = 0;
			continue;
		}
		//caso (0,1) (x,1)  
		if (((primo[i] == 'x') || (primo[i] == '0')) && (secondo[i] == '1')) {
			comandi[i] = 1;
			continue;
		}
		//caso (1,0) (1,x)  
		if ((primo[i] == '1') && ((secondo[i] == '0') || (secondo[i] == 'x'))) {
			comandi[i] = 2;
			continue;
		}
		//caso (1,1)  
		if ((primo[i] == '1') && (secondo[i] == '1')) {
			comandi[i] = 3;
			continue;
		}
	}

	return comandi;
}


vector<string> read_save_file(string nomefile) {

	/*
	ifstream t(nomefile);
	string str;

	t.seekg(0, std::ios::end);
	str.reserve(t.tellg());
	t.seekg(0, std::ios::beg);

	str.assign((std::istreambuf_iterator<char>(t)),
		std::istreambuf_iterator<char>());
	return str;*/

	ifstream file(nomefile);
	vector<string> sequences;
	int index = 0;
	//parsing entire file
	while (!file.eof()) {
		string temp;
		//if there's a '>' ignore all files until new line, where the actual string starts
		if (file.peek() == '>')
			file.ignore(numeric_limits<streamsize>::max(), '\r\n');
		for (;;) {
			if (file.peek() == '>')
				break;
			getline(file, temp);
			temp.erase(temp.length() - 2, 2);
			sequences[index] = sequences[index] + temp;
		}
		index++;
	}
	return sequences;
}