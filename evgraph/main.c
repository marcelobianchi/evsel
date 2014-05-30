#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include <interaction.h>
#include <continents.h>
#include <plates.h>

char *legends[] = { "X", "Y", "",                                   /* Generic XY-Graph */
					"Longitude", "Latitude", "World Map",           /* Map */
					"Longitude", "Depth (km)", "Longitude Section", /* LonSection */
					"Latitude", "Depth (km)", "Latitude Section",   /* LatSection */
					"Magnitude", "Count", "Magnitude Histogram"     /* Magnitude Histogram */
				};

typedef struct points {
	float *x;
	float *y;
	long n;
	int colx;
	int coly;
	float xmin;
	float xmax;
	float ymin;
	float ymax;
} POINTS;

typedef struct graphcontrol {
	int haspoints;
	int haslines;
	int hascontinents;
	int hasplates;
	int growperc;

	char **legend;

	/* Line Parameters */
	int fitline;
	float a;
	float b;

	float xmin;
	float xmax;
	float ymin;
	float ymax;
} GRAPHCONTROL;

/*
 * Processors
 */
float mulminus(float value) {
	return value * -1;
}

float takelog(float value) {
	return log10(value);
}

/*
 * POINTS struct methods
 */
void addPt(POINTS *p, float x, float y){
	p->x = realloc(p->x,sizeof(float *) * (p->n+1));
	p->y = realloc(p->y,sizeof(float *) * (p->n+1));
	p->x[p->n] = x;
	p->y[p->n] = y;
	if (p->n == 0) {
		p->xmax = p->xmin = x;
		p->ymin = p->ymax = y;
	} else {
		if (x > p->xmax) p->xmax = x;
		if (y > p->ymax) p->ymax = y;
		if (x < p->xmin) p->xmin = x;
		if (y < p->ymin) p->ymin = y;
	}
	p->n++;
}

POINTS *newP(int xcol, int ycol){
	POINTS *p = malloc(sizeof(POINTS));
	p->n = 0;
	p->colx = xcol;
	p->coly = ycol;
	p->x = NULL;
	p->y = NULL;
	p->xmax = 1.0;
	p->xmin = 0.0;
	p->ymin = 0.0;
	p->ymax = 1.0;
	return p;
}

void *cleanP(POINTS **p){
	if ((*p)->x) free((*p)->x);
	if ((*p)->y) free((*p)->y);
	(*p)->x = NULL;
	(*p)->y = NULL;
	(*p)->n = 0;
	free((*p));
	(*p) = NULL;
}

/*
 * Execution methods
 */
int parseLine(char *line, POINTS *p, float (*xprocessor)(float), float (*yprocessor)(float)) {
	char *fs = line;
	int currentcol = 0;

	int gotx = 0;
	int goty = 0;

	float x = 0.0;
	float y = 0.0;

	while (1) {
		char *rt;
		float value = strtod(fs, &rt);
		if (value == 0.0 && fs == rt) break;
		currentcol++;
		if (currentcol == p->colx) {
			x = value;
			gotx = 1;
		}
		if (currentcol == p->coly) {
			y = value;
			goty = 1;
		}
		if (gotx && goty) break;
		fs = rt;
	}

	if (!gotx || !goty) return 1;

	if (xprocessor != NULL) x = xprocessor(x);
	if (yprocessor != NULL) y = yprocessor(y);

	addPt(p,x,y);
	return 0;
}

void plot(GRAPHCONTROL *gr, POINTS *p) {
	cpgenv(gr->xmin, gr->xmax, gr->ymin, gr->ymax, 0, 0);

	cpglab(gr->legend[0], gr->legend[1], gr->legend[2]);

	/* Graphs */
	if (gr->haspoints) {
		cpgpt(p->n, p->x, p->y, 1);
	}

	if (strncmp(gr->legend[0], "Longitude", 9) == 0 ) {
		int i;

		if (gr->hascontinents) {
			cpgsci(2);
			cpgslw(5);
			for(i=0; i < ncontinentes; i++) {
				if (continentes[i][0] == -999 && continentes[i][1] == 999 ) {
					i++;
					cpgmove(continentes[i][0], continentes[i][1]);
					continue;
				}
				cpgdraw(continentes[i][0], continentes[i][1]);
			}
		}

		if (gr->hasplates) {
			cpgsci(3);
			cpgslw(3);
			for(i=0; i < nplates; i++) {
				if (plates[i][0] == -999 && plates[i][1] == 999 ) {
					i++;
					cpgmove(plates[i][0], plates[i][1]);
					continue;
				}
				if (fabs(plates[i][0] - plates[i-1][0]) > 180) {
					cpgmove(plates[i][0], plates[i][1]);
				}
				cpgdraw(plates[i][0], plates[i][1]);
			}
		}
		cpgsci(1);
		cpgslw(1);

	}

	if (gr->haslines) {
		cpgline(p->n, p->x, p->y);
	}

	return;
}

void showERR(char *message) {
	fprintf(stderr, message);
	return;
}

void rangeadjust(GRAPHCONTROL *gr, POINTS *p, float *x1, float *x2, float *y1, float *y2) {
	/* Axis */
	float xgrow = (p->xmax - p->xmin) * gr->growperc;
	float xmin = (p->xmin - xgrow);
	float xmax = (p->xmax + xgrow);
	float ygrow = (p->ymax - p->ymin) * gr->growperc;
	float ymin = (p->ymin - ygrow);
	float ymax = (p->ymax + ygrow);

	if ( (x1 != NULL) && (x2 != NULL) && (y1 != NULL) && (y2 != NULL) && (*x2 > *x1) && (*y2 > *y1)) {
		xgrow = (*x2 - *x1) * gr->growperc;
		xmin  = (*x1 - xgrow);
		xmax  = (*x2 + xgrow);
		ygrow = (*y2 - *y1) * gr->growperc;
		ymin  = (*y1 - ygrow);
		ymax  = (*y2 + ygrow);
	}

	if (xmin == xmax) {
		xmin -= 1;
		xmax += 1;
	}

	if (ymin == ymax) {
		ymin -= 1;
		ymax += 1;
	}

	gr->xmax = xmax;
	gr->xmin = xmin;
	gr->ymax = ymax;
	gr->ymin = ymin;
}

void order(float *a, float *b) {
	float c;
	if (*a > *b) {
		c = *b;
		*b = *a;
		*a = c;
	}
	return;
}

int main(int argc, char **argv) {
	GRAPHCONTROL gr = {
		0,   /* Has Points */
		0,   /* Has Lines */
		0,   /* Has Continents */
		0,   /* Has Plates */
		0.1, /* Grow % */
		&legends[0], /* Default legend set */
		0,   /* Fit line */
		0.0, /* a */
		0.0, /* b */
		0,   /* xmin */
		0,   /* xmax */
		0,   /* ymin */
		0    /* ymax */
	};

	int i, ID;
	POINTS *p = NULL;

	float (*xprocessor)(float) = NULL;
	float (*yprocessor)(float) = NULL;

	/*
	 * Command line arguments
	 */
	for (i = 1; i < argc; i++) {
		/* Major modes */
		if (strncmp(argv[i], "map", 3) == 0) {
			gr.legend = &legends[3];
			if (p != NULL) showERR("Cannot select MAP mode, only one mode is accepted.");
			p = newP(8,7);
		}
		if (strncmp(argv[i], "graph", 5) == 0) {
			gr.legend = &legends[0];
			if (p != NULL) showERR("Cannot select GRAPH mode, only one mode is accepted.");
			p = newP(1,2);
		}
		if (strncmp(argv[i], "lonsection", 10) == 0) {
			gr.legend = &legends[6];
			if (p != NULL) showERR("Cannot select lonsection mode, only one mode is accepted.");
			p = newP(8,9);
			yprocessor = &mulminus;
		}
		if (strncmp(argv[i], "latsection", 10) == 0) {
			gr.legend = &legends[9];
			if (p != NULL) showERR("Cannot select latsection mode, only one mode is accepted.");
			p = newP(7,9);
			yprocessor = &mulminus;
		}
		if (strncmp(argv[i], "maghist", 7) == 0) {
			gr.legend = &legends[12];
			if (p != NULL) showERR("Cannot select maghist mode, only one mode is accepted.");
			p = newP(10,-1);
		}

		/* Options */
		if (strncmp(argv[i], "dot", 3) == 0) {
			gr.haspoints = 1;
		}
		if (strncmp(argv[i], "lines", 5) == 0) {
			gr.haslines = 1;
		}
		if ((strncmp(argv[i], "continent", 9) == 0)||(strncmp(argv[i], "coast", 5) == 0)) {
			gr.hascontinents = 1;
		}
		if (strncmp(argv[i], "plates", 5) == 0) {
			gr.hasplates = 1;
		}
		if (strncmp(argv[i], "semilogy", 8) == 0) {
			yprocessor = &takelog;
		}
		if (strncmp(argv[i], "semilogx", 8) == 0) {
			xprocessor = &takelog;
		}
		if (strncmp(argv[i], "log", 3) == 0) {
			xprocessor = &takelog;
			yprocessor = &takelog;
		}
	}

	/* When no representation is selected haspoints is default to on */
	if ((gr.haslines == 0) && (gr.haspoints ==0)) {
		gr.haspoints = 1;
	}

	/*
	 * File Parsing
	 */
	char line[1024];
	fprintf(stderr,"Reading from STDIN ... ");
	while (fgets(line, 1024, stdin)) {
		parseLine(line, p, xprocessor, yprocessor);
	}
	fprintf(stderr, "%d points\n", p->n);

	/*
	 * Control Block
	 */
	ID = cpgopen("/xwindow");
	cpgask(0);

	char ch = 'Z';
	float x1, x2, y1, y2;
	rangeadjust(&gr, p, NULL, NULL, NULL, NULL);
	while (ch != 'Q') {
		plot(&gr, p);
		ch = getonechar(&x1, &y1, 0);
		switch (ch) {

		case('P'):
			cpgclos();
			cpgopen("/ps");
			plot(&gr, p);
			cpgclos();
			cpgopen("/xwindow");
			cpgask(0);
			break;

		case('O'):        /* Zoom */
			rangeadjust(&gr, p, NULL, NULL, NULL, NULL);
			break;

		case('A'):
		case('X'):
			x2 = x1;
			y2 = y1;
			getonechar(&x2, &y2, 1);
			order(&x1,&x2);
			order(&y1,&y2);
			rangeadjust(&gr, p, &x1, &x2, &y1, &y2);
			break;
		}
	}


	/*
	 * Clean up block
	 */
	cpgclos();

	cleanP(&p);
	if (p == NULL) fprintf(stderr,"Cleanned");

	/*
	 * Done
	 */
	return 0;
}
