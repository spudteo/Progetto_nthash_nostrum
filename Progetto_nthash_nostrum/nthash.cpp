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
static const uint64_t hash_NULLO = 0x0000000000000000;

vector<string> seqs;

uint64_t primo_Hash(string sequenza_input, int spaceSeedSize);
uint64_t restanti_Hash(uint64_t precedente, char uscente, char entrante, int spaceSeedSize);
vector <uint64_t> total_seq_hash(string sequenza_input, int spaceSeedSize);
vector<vector<uint64_t>> sequenzeTotaliSeedUni(string spacedSeed);
string complemento(string spacedSeed);
vector<int> preProcessing_1(string spacedSeedComp);
vector<int> preProcessing_2(string spacedSeedComp, vector<int> ab);
uint64_t toInt(char seq);
uint64_t leftRotate(uint64_t n, int d);
uint64_t rightRotate(uint64_t n, int d);
vector<string> read_save_file(string nomefile);
vector<string> getPastTok(int seqIndex, int tokIndex, vector<int> ab, int sslen);
uint64_t metodoNuovo_primiHash(string token, string spacedSeedComplementato);
uint64_t metodoNuovo_restantiHash(string past, string token, vector<int> ab, uint64_t hashVecchio, vector<int> comandi);
uint64_t hash_stupido(string sequenza_input, string spacedSeed);
vector<uint64_t> hash_stupido_UNASeq(string sequenza, string spacedSeed);
vector<vector<uint64_t>> hash_stupido_interaSeq(string spacedSeed);
vector<uint64_t> metodoNuovo_UNASeq(string sequenza, string spacedSeedComplementato, vector<int> ab, vector<int> comandi);
vector<vector<uint64_t>> metodoNuovo_interaSeq(string spacedSeedComplementato);
vector<vector<uint64_t>> finalHash(vector<vector<uint64_t>> complete, vector<vector<uint64_t>> spaced);


int main()
{
	//parse command line DA FARE//////////
	string nomefile = "C:\\Users\\teosp\\OneDrive\\Desktop\\teouni\\Algoritmi per la bioinformatica\\dataaset\\small_test.fna";
	//apre file 
	ifstream inputfile;
	inputfile.open(nomefile);
	seqs = read_save_file(nomefile);
	//non cancellare da qua in su se no non si puo fare un cazzo//


	string spacedSeed = "110010101101010101001";
	int spaceSeedSize = spacedSeed.length();

	//prova metodo "stupido"
	vector<vector<uint64_t>> prova_1 = hash_stupido_interaSeq(spacedSeed);

	vector<vector<uint64_t>> prova_2 = finalHash(sequenzeTotaliSeedUni(spacedSeed), metodoNuovo_interaSeq(spacedSeed));

	for (int i = 0; i < prova_1.size(); i++) {
		cout << " Hash sequenza: " << i << endl;
		for (int j = 0; j < prova_1.at(i).size(); j++) {
			cout << "Hash stupido: " << prova_1.at(i).at(j) << "  Hash nostro: " << prova_2.at(i).at(j) << endl;
		}
	}
	
	/*
	//prova dei primi hash con il nuovo metodo
	vector<string> token;
	token.push_back("ytgh");
	string sequenza;
	for (int i = 0; i < spacedSeed.length(); i++) {
		sequenza += seqs.at(0)[i];
	}
	token.push_back(sequenza);
	cout << "Hash nuovo metodo: " << metodoNuovo_primiHash(token,spacedSeed)<<endl;
	*/

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

	cout << "fine main" << endl;
	cin.ignore();
}

//SEED-tutti uni-trova il primo hash, dalla posizione 0 a k-1
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

//SEED-tutti uni-trova i restanti hash fino all'ultimo k-mer
uint64_t restanti_Hash(uint64_t precedente, char uscente, char entrante, int spaceSeedSize) {

	//hash di quello che sta prima ruotato di 1 XOR hash carattere che esce ruotato di k XOR hash di quello che esce 
	uint64_t hash = ((leftRotate(precedente, 1)) ^ (leftRotate(toInt(uscente), spaceSeedSize))) ^ (toInt(entrante));
	return hash;
}

//SEED-tutti uni-trova tutti gli hash e restituisce il vettore con gli hash, la posizione è l'indice di hash
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

//SEED-tutti uni-calcola tutti gli hash di tutte le sequenze che ci sono dentro seq e che legge dal file 
vector<vector<uint64_t>> sequenzeTotaliSeedUni(string spacedSeed) {
	vector<vector<uint64_t>> fullHash;
	for (int i = 0; i < seqs.size(); i++) {
		fullHash.push_back(total_seq_hash(seqs.at(i), spacedSeed.length()));
	}
	return fullHash;
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

//legge il file e mette le sequenze in un vettore di stringhe
vector<string> read_save_file(string nomefile) {

	ifstream file(nomefile);
	vector<string> sequenze;
	string temp;
	string temp2;
	int primo = 0;
	//legge una linea alla volta
	while (getline(file, temp)) {
		if (temp[0] == '>' || file.eof()) {
			if (primo > 0) {
				//inserisce la sequenza al vettore
				sequenze.push_back(temp2);
				temp2 = "";
			}
			primo++;
			continue;
		}
		else {
			//accumula le stringhe tra due caratteri >
			temp2.append(temp);
		}
	}
	return sequenze;
}

//restituisce il token e past 
vector<string> getPastTok(int seqIndex, int tokIndex, vector<int> ab, int sslen) {
	vector<string> pastTok (2);
	string seq = seqs[seqIndex];
	pastTok[1] = seq.substr(tokIndex, sslen);
	int dist = ab[0] - ab[1];
	if (tokIndex - dist < 0) {
		cout << "Past vector not available, returning empty past." << endl;
		pastTok[0] = "-1";
		return pastTok;
	}
	pastTok[0] = seq.substr(tokIndex - dist, sslen);
	return pastTok;
}

//calcola i primi (a-b) hash nel nostro nuovo metodo 
uint64_t metodoNuovo_primiHash(string token, string spacedSeedComplementato) {

	uint64_t hash = hash_NULLO;
	for (int i = 0; i < token.length(); i++) {
		if (spacedSeedComplementato[i] == '1') {
			hash = hash ^ leftRotate(toInt(token[i]), token.length() - 1 - i);
		}
	}
	return hash;
}

//calcola i restanti hash con il nostro nuovo metodo
uint64_t metodoNuovo_restantiHash(string past, string token, vector<int> ab, uint64_t hashVecchio,vector<int> comandi) {

	//rotate a sx di a 
	uint64_t newHash = leftRotate(hashVecchio,ab[0]);

	for (int i = 0; i < comandi.size(); i++) {
		switch (comandi[i])
		{
		//inserimento
		case 1:
			newHash = newHash ^ (leftRotate(toInt(token[i]), (comandi.size() - i - 1)));
			break;
		//cancellazione
		case 2:
			newHash = newHash ^ (leftRotate(toInt(past[i]), (comandi.size() - i - 1)));
			break;
		//entrambe
		case 3:
			newHash = newHash ^ (leftRotate(toInt(token[i]), (comandi.size() - i - 1)));
			newHash = newHash ^ (leftRotate(toInt(past[i]), (comandi.size() - i - 1)));
			break;
		//avevamo già il carattere nella posizione giusta
		default:
			break;
		}
	}

	//rotate a destra di b 
	newHash = rightRotate(newHash, ab[1]);

	return newHash;
}

//calcola tutti gli hash per una intera sequenza e li mette dentro un vettore che poi ritorna  
vector<uint64_t> metodoNuovo_UNASeq(string sequenza, string spacedSeedComplementato, vector<int> ab, vector<int> comandi) {

	vector<uint64_t> hashSeq;
	
	int primi = (ab[0] - ab[1]);
	//calcola i primi hash 
	for (int i = 0; i < primi; i++) {
		string token;
		for (int j = 0; j < spacedSeedComplementato.length(); j++) {
			token += sequenza[i + j];
		}
		hashSeq.push_back(metodoNuovo_primiHash(token,spacedSeedComplementato));
	}

	//calcoliamo i restanti hash
	//ogni possibile k sequenza dentro l'intera sequenza
	for (int i = primi; i <= sequenza.length() - spacedSeedComplementato.length(); i++) {
		//creo la stringa da dare in input al metodo che deve essere grande come lo spacedSeed
		string token;
		string past;

		for (int j = 0; j < spacedSeedComplementato.length(); j++) {
			token += sequenza[i + j];
		}
		for (int j = 0; j < spacedSeedComplementato.length(); j++) {
			past += sequenza[i-(primi) + j];
		}
		hashSeq.push_back(metodoNuovo_restantiHash(past, token, ab, hashSeq.at(i - (primi)), comandi));

	}


	return hashSeq; 
}

//calcola tutti gli hash di tutte le sequenze che ci sono dentro seq e che legge dal file
vector<vector<uint64_t>> metodoNuovo_interaSeq(string spacedSeed) {

	vector<vector<uint64_t>> fullHash;
	//preprocessing 
	string spacedSeedComplementato = complemento(spacedSeed);
	vector<int> ab = preProcessing_1(spacedSeedComplementato);
	vector<int> comandi = preProcessing_2(spacedSeedComplementato, ab);

	for (int i = 0; i < seqs.size(); i++) {
		fullHash.push_back(metodoNuovo_UNASeq(seqs.at(i), spacedSeedComplementato,ab,comandi));
	}
	return fullHash;
}

//calcola l'hash in maniera "stupida" andando a vedere ogni 0 ed ogni 1 nello spaced seed e fa i conti di conseguenza
uint64_t hash_stupido(string sequenza_input, string spacedSeed) {

	uint64_t hash = hash_NULLO;
	for (int i = 0; i < spacedSeed.length(); i++) {
		if (spacedSeed[i] == '1') {
			hash = hash ^ leftRotate(toInt(sequenza_input[i]), sequenza_input.length() - 1 - i);
		}
	}
	return hash;
}

//calcola tutti gli hash per una intera sequenza e li mette dentro un vettore che poi ritorna 
vector<uint64_t> hash_stupido_UNASeq(string sequenza, string spacedSeed) {

	vector<uint64_t> hashSeq;
	for (int i = 0; i <= sequenza.length() - spacedSeed.length(); i++) {
		//creo la stringa da dare in input al metodo che deve essere grande come lo spacedSeed
		string sequenza_input;
		for (int j = 0; j < spacedSeed.length(); j++) {
			sequenza_input += sequenza[i+j];
		}
		hashSeq.push_back(hash_stupido(sequenza_input,spacedSeed));
	}

	return hashSeq;
}

//calcola tutti gli hash di tutte le sequenze che ci sono dentro seq e che legge dal file 
vector<vector<uint64_t>> hash_stupido_interaSeq(string spacedSeed) {

	vector<vector<uint64_t>> fullHash;
	for (int i = 0; i < seqs.size(); i++) {
		fullHash.push_back(hash_stupido_UNASeq(seqs.at(i),spacedSeed));
	}
	return fullHash;
}

//calcola la differenza tra l'hash completo e quello ottenuto con lo spaced seeds
vector<vector<uint64_t>> finalHash(vector<vector<uint64_t>> complete, vector<vector<uint64_t>> spaced) {

	vector<vector<uint64_t>> finalH;

	for (int i = 0; i < complete.size(); i++) {
		vector<uint64_t> uniMenoSpace;
		for (int j = 0; j < complete.at(i).size(); j++) {
			uniMenoSpace.push_back(complete.at(i).at(j) ^ spaced.at(i).at(j));
		}
		finalH.push_back(uniMenoSpace);
	}
	return finalH;
}
