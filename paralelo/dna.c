#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <omp.h>
#include <time.h>

// MAX char table (ASCII)
#define MAX 256

// Boyers-Moore-Hospool-Sunday algorithm for string matching
int bmhs(char *string, int n, char *substr, int m) {

	int d[MAX];
	int i, j, k;
	// pre-processing 
	#pragma omp simd //tacerto
	for (j = 0; j < MAX; j++) // talvez de pra usar o pragm do smid aqui (so no primeiro, o segundo da pq os elementos nao estao anihados)
		d[j] = m + 1;
	
	#pragma omp simd
	for (j = 0; j < m; j++)
		d[(int) substr[j]] = m - j;
	// searching
	i = m - 1;
	while (i < n) {
		k = i;
		j = m - 1;
		while ((j >= 0) && (string[k] == substr[j])) {
			j--;
			k--;
		}
		if (j < 0)	
		{
			return k + 1;
		}
		i = i + d[(int) string[i + 1]];
	}
	// printf("");
	return -1;
}

FILE *fdatabase, *fquery, *fout;

//[tempo e curto pra krl, nao tem necessidade de usar isso] fazer aquele esquema do pragma de single thread e transformar cada if em uma task... talvez pro close files tb seja uma boa ideia (teria que medir o tempo)
void openfiles() {

	fdatabase = fopen("dna.in", "r+");
	if (fdatabase == NULL) {
		perror("dna.in");
		exit(EXIT_FAILURE);
	}

	fquery = fopen("query.in", "r");
	if (fquery == NULL) {
		perror("query.in");
		exit(EXIT_FAILURE);
	}

	fout = fopen("dna.out", "w");
	if (fout == NULL) {
		perror("fout");
		exit(EXIT_FAILURE);
	}

}

void closefiles() {
	fflush(fdatabase);
	fclose(fdatabase);

	fflush(fquery);
	fclose(fquery);

	fflush(fout);
	fclose(fout);
}

static inline void remove_eol(char *line) {

	int i = strlen(line) - 1;
	while (line[i] == '\n' || line[i] == '\r') {
		line[i] = 0;
		i--;
	}
}

typedef struct{
	char text_cab[100];
	char genoma[1000001];
} marcador;

marcador marcadores[1000]; // 1000 marcadores inicialmente, acho que nao passa disso
int cabs;


void preprocessamento(){
//pegar a primeira linha ( cabecalho inicial) fora do while, e marcar a pulada de linha dele com um &
//depois trabalhar em cima desse simbolo
// os marcadores [i][i+1] sao representantes das posicoes de real inicio/parada para cada thread
// na hora de quebrar os lacos, cada thread via ter que ser seu numero
// pq nao esta copiando depois do primeiro /n?
// se for assim com um campo so para texto, entao vai mudar os marcadores. ai vai mudar os controladores de id
// ta errando so a onde pegar o paraneu
	char linha[100];
	int j; //i;
	fseek(fdatabase, 0, SEEK_SET);
	fgets(linha, 100, fdatabase);
	remove_eol(linha);

	strcpy( marcadores[0].text_cab, linha );
	j = 0; //i = 0;

	while(!feof(fdatabase))
	{
		fgets(linha, 100, fdatabase);
		remove_eol(linha);
		do { 
				strcat(marcadores[j].genoma, linha);
				
				if (fgets(linha, 100, fdatabase) == NULL)
					break;
				remove_eol(linha);
		} while (linha[0] != '>');
		cabs++;
		j++;
		strcpy( marcadores[j].text_cab, linha );
	}
	cabs-=1;
}

char *bases;
char *str;

int main(char** argv) {
	clock_t t = clock();
	cabs =0;
	int nt = 0;// (int) *argv[1]; // ver outra forma de usar isso 
	bases = (char*) malloc(sizeof(char) * 1000001);
	if (bases == NULL) {
		perror("malloc");
		exit(EXIT_FAILURE);
	}
	str = (char*) malloc(sizeof(char) * 1000001);
	if (str == NULL) {
		perror("malloc str");
		exit(EXIT_FAILURE);
	}

	openfiles();

	// printf("bacs antes do openfile: %d \n",cabs);
	preprocessamento();
	// printf("bacs dep do openfile: %d \n",cabs);
	// for (int m = 0; ( marcadores[m].cab != -1 || marcadores[m].com != -1 ); m++ )
	// 	printf("posicao e valor: [ patogenico: %s \n ] [ m: %d ] ~ [cab: %d] ~ [com %d] ~ [elemento da pos: %c  %c] \n", marcadores[m].text_cab, m, marcadores[m].cab, marcadores[m].com, dna_completo[ (int)marcadores[m].cab ], dna_completo[ (int) marcadores[m].com] );
	// for (int y = 0; y < cabs; y++)
	// {	
	// 	printf("seq: %s \n tenho tamanho: %ld \n %.30s \n", marcadores[y].text_cab, strlen( marcadores[y].genoma ) ,marcadores[y].genoma);
	// }
	// return 0;
	char desc_dna[100], desc_query[100];
	char line[100];
	int i, found, result;
	int y,h=0;
	fgets(desc_query, 100, fquery); 
	remove_eol(desc_query);
	while (!feof(fquery)) {
		fprintf(fout, "%s\n", desc_query);
		// read query string
		fgets(line, 100, fquery); // vai precisar de barreira 
		remove_eol(line);
		str[0] = 0;
		i = 0;
		do {  // talvez mudar esse bloco para algo que seja possivel fazer de forma paralela, 2 linhas por vez ou algo do tipo
			strcat(str + i, line); // vai precisar de barreira
			if (fgets(line, 100, fquery) == NULL)
				break;
			remove_eol(line);
			i += 80;
		} while (line[0] != '>');
		strcpy(desc_query, line); // vai precisar de barreira

		// read database and search
		found = 0;
		fseek(fdatabase, 0, SEEK_SET);
		fgets(line, 100, fdatabase);
		remove_eol(line);
		
		// printf("antes do for: %d \n", h);
		for (y = 0; y <=cabs; y++ )
		{	
			// printf("entrei no for meu valor> %d ||| hzao: %d \n", y, h );

			// printf("pos result \n");
			result = bmhs( marcadores[y].genoma, strlen( marcadores[y].genoma ), str, strlen(str)); //str eh a query
			// printf(" ** entrei pos for meu valor> %d ||| hzao: %d \n", y, h );
			if (result > 0) {
				// printf("antes de escrever o result");
				fprintf(fout, "%s\n%d\n", marcadores[y].text_cab, result);
				found++;
				// printf("depois de escrever o result");
			}
			// printf("pos if\n");
		}
		if (!found)
			fprintf(fout, "NOT FOUND\n");
		h++;
	}

	closefiles();

	free(str);
	free(bases);
	t = clock() - t; //dado em seg
	double t_total = ((double)t)/CLOCKS_PER_SEC;
	printf("Resultado openfile: %f \n", t_total);
	return EXIT_SUCCESS;
}





	// t = clock() - t; //dado em seg
	// double t_total = ((double)t1)/CLOCKS_PER_SEC;
	// printf("Resultado openfile: %f \n", t_total);