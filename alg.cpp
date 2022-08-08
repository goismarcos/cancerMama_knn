#include <iostream>
#include <vector>
#include <math.h>
#include <set>
#include <map>

using namespace std;

//classe  célula
class Celula
{
	private:
		/*
		PARAMETROS USADOS NO DATASET
		1. Número do código da amostra: número de identificação
		2. Clump Espessura: 1 - 10
		3. Uniformidade do tamanho da célula: 1 - 10
		4. Uniformidade da Forma da Célula: 1 - 10
		5. Adesão Marginal: 1 a 10
		6. Tamanho Único de Células Epiteliais: 1 - 10
		7. Núcleos Nus: 1 - 10
		8. Cromatina Branda: 1 a 10
		9. Nucleoli Normal: 1 - 10
		10. Mitoses: 1 - 10
		11. Classe: (2 para benigno, 4 para maligno)
		*/
		//parametros célula
		double  clump_espessura,
				uniformidade_tamanho_celula,
				uniformidade_forma_celula,
				adesao_marginal,
				tamanho_unico_celulas_epiteliais,
				nucleos,
				cromatina_branda,
				nucleoli_normal,
				mitoses;

		string	classe;
	public:
		//método construtor para receber os parametro
		Celula(double clump_espessura, double uniformidade_tamanho_celula,
			   double uniformidade_forma_celula, double adesao_marginal,
			   double tamanho_unico_celulas_epiteliais, double nucleos,
			   double cromatina_branda, double nucleoli_normal,
			   double mitoses, string classe)
			   {
			   		this->clump_espessura = clump_espessura;
			   		this->uniformidade_tamanho_celula = uniformidade_tamanho_celula;
			   		this->uniformidade_forma_celula = uniformidade_forma_celula;
			   		this->adesao_marginal = adesao_marginal;
			   		this->tamanho_unico_celulas_epiteliais = tamanho_unico_celulas_epiteliais;
			   		this->nucleos = nucleos;
			   		this->cromatina_branda = cromatina_branda;
			   		this->nucleoli_normal = nucleoli_normal;
			   		this->mitoses = mitoses;
			   		this->classe = classe;
			   }

		//metodos get para obter parametros
		double get_clump_espessura()
		{
			return clump_espessura;
		}

		double get_uniformidade_tamanho_celula()
		{
			return uniformidade_tamanho_celula;
		}

		double get_uniformidade_forma_celula()
		{
			return uniformidade_forma_celula;
		}

		double get_adesao_marginal()
		{
			return adesao_marginal;
		}

		double get_tamanho_unico_celulas_epiteliais()
		{
			return tamanho_unico_celulas_epiteliais;
		}

		double get_nucleos()
		{
			return nucleos;
		}

		double get_cromatina_branda()
		{
			return cromatina_branda;
		}

		double get_nucleoli_normal()
		{
			return nucleoli_normal;
		}

		double get_mitoses()
		{
			return mitoses;
		}

		string get_classe()
		{
			return classe;
		}

};//fim da classe celula

//função para calcular a distancia eu clidiana
double Distancia_euclidiana(Celula celula_1, Celula celula_2 )
{

	double resultado_soma =  pow((celula_1.get_clump_espessura() - celula_2.get_clump_espessura()), 2) +
							 pow((celula_1.get_uniformidade_tamanho_celula() - celula_2.get_uniformidade_tamanho_celula()), 2) +
							 pow((celula_1.get_uniformidade_forma_celula() - celula_2.get_uniformidade_forma_celula()), 2) +
							 pow((celula_1.get_adesao_marginal() - celula_2.get_adesao_marginal()), 2) +
							 pow((celula_1.get_tamanho_unico_celulas_epiteliais() - celula_2.get_tamanho_unico_celulas_epiteliais()), 2) +
							 pow((celula_1.get_nucleos() - celula_2.get_nucleos()), 2) +
							 pow((celula_1.get_cromatina_branda() - celula_2.get_cromatina_branda()), 2) +
							 pow((celula_1.get_nucleoli_normal() - celula_2.get_nucleoli_normal()), 2) +
							 pow((celula_1.get_mitoses() - celula_2.get_mitoses()), 2) ;

	//retornando a raiz quadrada da soma da diferença dos quadrados
	return sqrt(resultado_soma);
}

//funcao de classificação de amostras
//Esta função recebera 3 parametros o primeiro sera um vetor de celulas, o segundo sera um novo exemplo de celula,
// e o terceiros sera uma k que determinara quantos vizinhos mais proximos eu vou pegar.
string Classificador_amostras_celulas(vector<Celula>& celulas, Celula novo_exemplo_celula, int k)
{
	//não é recomendado que k seje par, pois pode ocorrer de ter o mesmo numero de classe para a celula
	//entao forçarei o k a ser impar decrementando ele, e maior que 0
	if( k % 2 == 0)
	{
		k--;
		if(k <= 0)
			k = 1;
	}

	int tamanho_vetor_de_celulas = celulas.size();

	//set de pear da distancia de cada celula do conjunto de treinamento para o novo exemplo de celula
	//cada pear vai ser composto pela distancia e pelo indice do vetor
	set<pair <double, int> > distancia_celulas;

	//calculando a distancia euclidiana do novo exemplo da celula, para cada amostra do conjunto de treinamento
	for(int i = 0; i < tamanho_vetor_de_celulas; i++)
	{
		double distancia = Distancia_euclidiana(celulas[i], novo_exemplo_celula);

		//enserindo no pear a distancia e o indice do vetor
		distancia_celulas.insert(make_pair(distancia,i));
	}

	//para decidir a qual classe pertence o novo exemplo de celula, basta verificar a classe mais frenquente entre os ks vizinhos mais proximos
	//do conjunto de treinamento.
	//para isso vou precisar de um set de pear iterator, e tambem um vetor para cada classe fazer a contagem de celula benigna ou maligna

	set<pair<double, int> >::iterator it;

	//declarando um vetor para as duas classes, posição 0 = benigna, posição 1 = maligna
	vector<int> contador_classes(2);

	int contador_k = 0;

	//percorrendo o set de pear das distancias das celulas
	for(it = distancia_celulas.begin(); it != distancia_celulas.end(); it++)
	{
		contador_k++;//incrementando o contador para chegar nos ks visinhos

		//ordem de parada quando o contador for igual ao numero de ks visiznhos
		if(contador_k == k)
			break;

		//pegando a classe da celula
		string classe = celulas[it->second].get_classe();

		//incrementando o contador da classe
		if(classe == "2")
			contador_classes[0]++;
		else if(classe == "4")
			contador_classes[1]++;
	}

	//agora vamos descobrir qual a classe mais frequente para retornar ela
	string classe_classificacao;

	//se a classe benigna for maior que a maligna, no primeiro if eu coloco igual, mais nunca ira acontecer de
	//ocorrer empate de classes pois meu k sempre sera impar e so temos duas classes
	if(contador_classes[0] >= contador_classes[1])
		classe_classificacao = "2";
	else
		classe_classificacao = "4";

	//retornando a classe predominante
	return classe_classificacao;
}

int main(int argc, char *argv[])
{
	setlocale(LC_ALL, "Portuguese");
	//declarando um vetor de celulas
	vector<Celula> celulas;

	//k que determinara a quantidade de vizinhos proximos
	int k = 5;

	//variavel para mostrar a precisao do teste
	float precisao = 0;

	//variavel do tamanho do conjunto de dados de treinamento
	//normalmente eé colocado 75 % do conjunto de dados para treinamento e os outros 25% para testes
	int tamanho_treinamento = 511;

	//o processo de treinamento consiste em armazenar o conjunto de dados de treinamento
	//e é o que vamos fazer agora
	for(int i = 0; i < tamanho_treinamento; i++)
	{
		//declarando as variaveis
		double  clump_espessura,
				uniformidade_tamanho_celula,
				uniformidade_forma_celula,
				adesao_marginal,
				tamanho_unico_celulas_epiteliais,
				nucleos,
				cromatina_branda,
				nucleoli_normal,
				mitoses;

		string classe;

		//fazendo a entrada de dados
		cin >> clump_espessura
			>> uniformidade_tamanho_celula
			>> uniformidade_forma_celula
			>> adesao_marginal
			>> tamanho_unico_celulas_epiteliais
			>> nucleos
			>> cromatina_branda
			>> nucleoli_normal
			>> mitoses
			>> classe;

			//inserindo no do vetor os paramentros de entrada
			celulas.push_back(Celula(	clump_espessura,
										uniformidade_tamanho_celula,
										uniformidade_forma_celula,
										adesao_marginal,
										tamanho_unico_celulas_epiteliais,
										nucleos,
										cromatina_branda,
										nucleoli_normal,
										mitoses,
										classe
									)
							 );
	}

	//numeros de acertos
	int acertos = 0;

	//variavel do tamanho do conjunto de dados de teste 25% , 683 representa o total de amostras que tenho no dataset
	int tamanho_teste = 683 - tamanho_treinamento ;

	//muitos importante salientar que esses dados de teste nao foram usados nos dados de treinamento, eles sao exlusivos
	//como se nunca foram usados

	//aogora fazeremos o processo de classificação desses teste
	for(int i = 0; i < tamanho_teste; i++)
	{
		//declarando as variaveis
		double  clump_espessura,
				uniformidade_tamanho_celula,
				uniformidade_forma_celula,
				adesao_marginal,
				tamanho_unico_celulas_epiteliais,
				nucleos,
				cromatina_branda,
				nucleoli_normal,
				mitoses;
		string classe;

		//fazendo a entrada de dados
		cin >> clump_espessura
			>> uniformidade_tamanho_celula
			>> uniformidade_forma_celula
			>> adesao_marginal
			>> tamanho_unico_celulas_epiteliais
			>> nucleos
			>> cromatina_branda
			>> nucleoli_normal
			>> mitoses
			>> classe;

		//classe celula_teste recebendo os parametros de entrada
		Celula celula_teste(clump_espessura,
							uniformidade_tamanho_celula,
							uniformidade_forma_celula,
							adesao_marginal,
							tamanho_unico_celulas_epiteliais,
							nucleos,
							cromatina_branda,
							nucleoli_normal,
							mitoses,
							classe
						   );

		string classe_obtida = Classificador_amostras_celulas(celulas, celula_teste, k);

		//mostrando classe espera e classe obtida
		//cout << i <<"° - " << classe <<" -> " << classe_obtida << "\n";

		//fazendo a contagem de acertos
		if(classe == classe_obtida)
		{
			cout << i <<"° - " << classe <<" -> " << classe_obtida << " *\n";
			acertos++;
		}
		else
		    cout << i <<"° - " << classe <<" -> " << classe_obtida << " x\n";



	}

	cout <<"----------------------------------------------------------------------------";
	//mostrando total de acertos e o total de testes feitos
	cout << "\n" <<acertos << " acertos de um total de " << tamanho_teste << " testes!\n";

	//calculando precisao do teste
	float aux;
	aux = acertos;
	precisao = ((aux/tamanho_teste)*100);

	//mostrando a porcentagem da precisao dos testes feitos
	cout << "A precisão de acertos foi de : "<< precisao << "%\n";
	cout <<"----------------------------------------------------------------------------\n";

	return 0;
}

