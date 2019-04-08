#include <iostream>
#include <string>
#include <vector>

using namespace std;
#define LONG_LONG_BITS 64

// 64-bit random seeds corresponding to bases and their complements
static const uint64_t seedA = 0x3c8bfbb395c60474;
static const uint64_t seedC = 0x3193c18562a02b4c;
static const uint64_t seedG = 0x20323ed082572324;
static const uint64_t seedT = 0x295549f54be24456;

uint64_t primo_Hash(string sequenza_input,int spaceSeedSize);
uint64_t restanti_Hash(uint64_t precedente, char uscente, char entrante, int spaceSeedSize );
vector <uint64_t> total_seq_hash(string sequenza_input, int spaceSeedSize);
uint64_t hash_complemento(string spacedSeed, string sequenza_input);
int preProcessing(string spacedSeed);
uint64_t toInt(char seq);
uint64_t leftRotate(uint64_t n, int d);
string complementary(string spacedSeed);

//fare una prova dei metodi fatti
int main()
{
	string spacedSeed = "1010";
	int spaceSeedSize = spacedSeed.length();
	
	vector<uint64_t> hashVector;
	string prova = "ACGTACGTACTGTCAGTCGATGCCTAGGTCATGCATGCGTCATGCATACTGACTATCTACTACTACTATCACC";
	hashVector = total_seq_hash(prova,spaceSeedSize);
	for (int i = 0; i < hashVector.size(); i++) {
		cout << hashVector[i] << endl;
	}

	hash_complemento(spacedSeed,prova);

	cout << complementary(spacedSeed) << endl;

	cout << preProcessing(complementary(spacedSeed)) << endl;

	cin.ignore();
}

//trova il primo hash, dalla posizione 0 a k-1
uint64_t primo_Hash(string sequenza_input,int spaceSeedSize) {

	//il primo carattere
	uint64_t hash= leftRotate(toInt(sequenza_input[0]), spaceSeedSize - 1);
	uint64_t temp;
	//ciclo lungo come il k-mer
	for (int i = 1; i < spaceSeedSize; i++)
	{
		//calcola lo shift necessario e poi fa lo XOR con quanto già calcolato 
		temp = leftRotate(toInt(sequenza_input[i]),spaceSeedSize -i +1);
		hash = hash ^ temp;
	}
	
	return hash;
} 

//trova i restanti hash fino all'ultimo k-mer
uint64_t restanti_Hash(uint64_t precedente,char uscente, char entrante,int spaceSeedSize) {

	//hash di quello che sta prima ruotato di 1 XOR hash carattere che esce ruotato di k XOR hash di quello che esce 
	uint64_t hash = leftRotate(precedente, 1) ^ leftRotate(toInt(uscente), spaceSeedSize) ^ toInt(entrante);
	return hash;
}

//trova tutti gli hash e restituisce il vettore con gli hash, la posizione è l'indice di hash
vector <uint64_t> total_seq_hash(string sequenza_input, int spaceSeedSize) {

	//dentro il vettore ci salvo tutti gli hash 
	vector <uint64_t> hashVector (sequenza_input.length() - spaceSeedSize +1);
	//per fare il primo hash mi servono solo k valori, taglio la sequenza di conseguenza
	string primicaratteri = sequenza_input.substr(0,spaceSeedSize-1);
	//calcolo l'hash del primo elemento e lo salvo
	uint64_t primoHash = primo_Hash(primicaratteri, spaceSeedSize);
	hashVector[0] = primoHash;

	//calcolo tutti gli altri hash e li aggiungo nel vettore
	//indice parte da 1 perchè 0 già calcolato,
	//si ferma ad i = (l-k) (lunghezza input- lunghezza finestra)
	for (int i = 1; i <= sequenza_input.length() - spaceSeedSize ; i++)
	{
		//hash di quello che sta prima, il primo carattere che esce, il primo carattere che entra nella finestra lunga k  
		hashVector[i] = restanti_Hash(hashVector[i-1],sequenza_input[i], sequenza_input[spaceSeedSize +i -1],spaceSeedSize);
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
	/* In n<<d, last d bits are 0. To
	put first 3 bits of n at
	last, do bitwise or of n<<d
	with n >>(INT_BITS - d) */
	return (n << d) | (n >> (LONG_LONG_BITS - d));
}



///////////////////////CAAPIRE COME FARE HASH DEL SINGOLO E BASTA//////////forse non serve///
//calcola la stringa complemento della stringa spaceSeed
uint64_t hash_complemento(string spacedSeed,string sequenza_input)
{
	//calcolo il complemento che mi serve 
	int pos1=0;
	int pos2=0;
	bool primo = true;
	for (int i = 0; i < spacedSeed.length(); i++) {
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
	for (int i = pos1; i <=pos2; i++) {
		if (spacedSeed[i] == '0') 
			complemento.append("1");
		if (spacedSeed[i] == '1') 
			complemento.append("0");
	}
	
	//calcolo l'hash della stringa complemento
	uint64_t hash_complemento = 4445;
	return hash_complemento;
}

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

////preprocessing/////
int preProcessing(string spacedSeedComp) {

	//posizione dove conviene attaccarsi
	int pos = 0;
	int score = 0;
	int bestScore = score;
	string temp = spacedSeedComp;
	//toglie ultimo carattere e inserisce uno 0 in prima posizione
	temp.pop_back();
	temp.insert(0, 1, '0');
	int iter = 1;
	while (iter < spacedSeedComp.length()) {
		//calcolo lo score di questo "attacco" di spacedseed
		for (int j = 0; j < spacedSeedComp.length(); j++) {
			if ((temp[j] == '1') && (spacedSeedComp[j] == '1'))
				score++;
		}
		if (score > bestScore) {
			bestScore = score;
			pos = iter;
		}
		score = 0;
		//riaggiusto la stringa 
		temp.pop_back();
		temp.insert(0, 1, '0');
		iter++;
	}
	return pos;

}

reffrfeferfe