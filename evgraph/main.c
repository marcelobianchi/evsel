#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <interaction.h>
#include <continents.h>
#include <plates.h>
#include <math.h>

char *legends[] = { "X", "Y", "",                      /* Generic XY-Graph */
                    "Longitude", "Latitude", "Map",    /* Map */
                    "Inline", "Depth (km)", "Section", /* Section */
                    "Mag", "Count", "Histogram"        /* Magnitude Histogram */
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

void addPt(POINTS *p, float x, float y){
    p->x = realloc(p->x,sizeof(float) * (p->n+1));
    p->y = realloc(p->y,sizeof(float) * (p->n+1));
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

int parseLine(char *line, POINTS *p) {
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

    addPt(p,x,y);
    return 0;
}

POINTS *newP(){
    POINTS *p = malloc(sizeof(POINTS));
    p->n = 0;
    p->colx = 1;
    p->coly = 2;
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

void plot(GRAPHCONTROL *gr, POINTS *p) {
    if (strncmp(gr->legend[1], "Depth", 5) == 0 ) {
        cpgenv(gr->xmin, gr->xmax, gr->ymax, gr->ymin, 0, 0);
    } else {
        cpgenv(gr->xmin, gr->xmax, gr->ymin, gr->ymax, 0, 0);
    }

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

int main(int argc, char **argv) {
    GRAPHCONTROL gr = {
        1,
        0,
        0,
        0,
        0.2,
        &legends[0],
        0,
        0.0,
        0.0, 0, 0, 0, 0
    };

    POINTS *p = newP();
    int i;

    int ID = cpgopen("/xwindow");
    cpgask(0);

    /* Command line arguments */
    for (i = 1; i < argc; i++) {
        if (strncmp(argv[i], "map", 3) == 0) {
            gr.legend = &legends[3];
        }
        if (strncmp(argv[i], "graph", 5) == 0) {
            gr.legend = &legends[0];
        }
        if (strncmp(argv[i], "section", 7) == 0) {
            gr.legend = &legends[6];
        }
        if (strncmp(argv[i], "maghist", 7) == 0) {
            gr.legend = &legends[9];
        }
        if (strncmp(argv[i], "lines", 5) == 0) {
            gr.haslines = 1;
        }
        if (strncmp(argv[i], "continents", 5) == 0) {
            gr.hascontinents = 1;
        }
        if (strncmp(argv[i], "plates", 5) == 0) {
            gr.hasplates = 1;
        }
    }

  char line[1024];
  fprintf(stderr,"Reading from STDIN ... \n");
  while (fgets(line, 1024, stdin)) {
     parseLine(line, p);
  }

//    for(i=0;i<34;i++) {
//        float r = random();
//        addPt(p, i, i*i + (r / RAND_MAX) * 45.0);
//    }

    char ch = 'Z';
    float x1,x2, y1,y2;
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
            rangeadjust(&gr, p, &x1, &x2, &y1, &y2);
            break;
        }
    }

    fprintf(stderr, "Read %d points", p->n);
    cleanP(&p);

    if (p == NULL) fprintf(stderr,"Cleanned");

    cpgclos();

    return 0;
}
