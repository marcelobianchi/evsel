/*
    This file is part of evsel package.
    
    evgraph is a tool to map 2D point data as extracted by the evsel
    SHELL script package.

    Copyright (C) 2013  Marcelo B. de Bianchi

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include <interaction.h>
#include <continents.h>
#include <plates.h>
#include <X11/Xlib.h>

#define FS 0.7

#define COLORNONE 0
#define COLORDEPTH 1
#define COLORMAG 2

#define ON 3
#define OFF 1

#define LAT 0
#define LON 1
#define DEP 2
#define MAG 3

/*
 * Structures
 */

typedef struct set {
	float *x;
	float *y;
	float *m;
	float *d;

	int y1;
	int y2;
	float d1;
	float d2;
	float m1;
	float m2;
	float lat1;
	float lat2;
	float lon1;
	float lon2;

	int region;

	long n;
} SET;

typedef struct graphcontrol {
	int haspoints;
	int haslines;
	int hascontinents;
	int hasplates;
	int growperc;
	int colormode;

	/* Line Parameters */
	int fitline;
	float a;
	float b;

	float xmin;
	float xmax;
	float ymin;
	float ymax;
} GRAPHCONTROL;

SET mainset;

GRAPHCONTROL control = {
	1,   /* Has Points */
	0,   /* Has Lines */
	1,   /* Has Continents */
	1,   /* Has Plates */
	0.1, /* Grow % */
	COLORDEPTH,
	0,   /* Fit line */
	0.0, /* a */
	0.0, /* b */
	0,   /* xmin */
	0,   /* xmax */
	0,   /* ymin */
	0    /* ymax */
};

/*
 * Color derivation function
 */

int magcolor(float mag) {
	if (mag < 4.)
		return 1;
	else if (mag< 6)
		return 2;
	else
		return 3;
}

int depthcolor(float depth) {
	if (depth < 35.)
		return 1;
	else if (depth < 80)
		return 2;
	else if (depth < 120)
		return 3;
	else if (depth < 300)
		return 4;
	else
		return 5;
}

/*
 *Graphic control & aux functions
 */

void resizemax(float scale)
{
	Display *disp;
	float ax, ay;
	int X, Y;

	/* Xlib code */
	disp = XOpenDisplay(NULL);
	if (disp == NULL) {
		fprintf(stderr, "No Display.\n");
		exit(-1);
	} else {
		Y = XDisplayHeightMM(disp, 0);
		X = XDisplayWidthMM(disp, 0);
	}
	XCloseDisplay(disp);
	/* End of Xlib code */

	ay = (double) Y / (double) X;
	ax = X / 25.4 * scale;
	cpgpap(ax, ay);
	cpgpage();
}

void minmax(float *x, int n, float *min, float *max) {
	int i;
	*min = *max = x[0];
	for (i = 0; i < n; i++) {
		if (*min > x[i]) *min = x[i];
		if (*max < x[i]) *max = x[i];
	}
	return;
}

void rangeadjust(GRAPHCONTROL *gr, SET *p, float *x1, float *x2, float *y1, float *y2) {
	/* Axis */
	float xgrow = (p->lon2 - p->lon1) * gr->growperc;
	float xmin = (p->lon1 - xgrow);
	float xmax = (p->lon2 + xgrow);
	float ygrow = (p->lat2 - p->lat1) * gr->growperc;
	float ymin = (p->lat1 - ygrow);
	float ymax = (p->lat2 + ygrow);

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

void orderint(int *a, int *b) {
	int c;
	if (*a > *b) {
		c = *b;
		*b = *a;
		*a = c;
	}
	return;
}

float linefit(float *x,float *y,int npts,float *a,float *b)
{
	int i;
	float D,Sxx,Sy,Sx,Sxy;
	float rms=0.0;

	Sxx=Sx=Sy=Sxy=D=0;
	for(i=0;i<npts;i++)
	{
		Sxx+=powf(x[i],2);
		Sx+=x[i];
		Sy+=y[i];
		Sxy+=x[i]*y[i];
	}
	D=npts*Sxx-powf(Sx,2);
	*b=(1.0/D)*(Sxx*Sy-Sx*Sxy);
	*a=(1.0/D)*((float)npts*Sxy-Sx*Sy);

	for(i=0;i<npts;i++)
		rms+=pow(y[i]-((*a)*x[i]+(*b)),2);

	rms=sqrt(rms/npts);
	return rms;
}

/*
 * Data Structure handling
 */

void load(SET *s) {
	char cmd[1024];

//	fprintf(stderr,"Start load.");
	if (s->n!=0) {
		if (s->x != NULL) free(s->x);
		if (s->y != NULL) free(s->y);
		if (s->m != NULL) free(s->m);
		if (s->d != NULL) free(s->d);

		s->x = NULL;
		s->y = NULL;
		s->m = NULL;
		s->d = NULL;

		s->n = 0;
	}

//	fprintf(stderr,"CMD LINE.");
	sprintf(cmd,"./run.sh  %d %d %f %f %f %f %f %f %f %f _TMP",
			s->y1,
			s->y2,
			s->m1,
			s->m2,
			s->d1,
			s->d2,
			s->lat1,
			s->lat2,
			s->lon1,
			s->lon2);
	system(cmd);

//	fprintf(stderr,"PROPER LOAD.");
	FILE *ent;
	ent = fopen("_TMP","r");
	if (ent == NULL) return;
	int i = 0;
	while (!feof(ent)) {
		s->x = (float*)realloc(s->x, sizeof(float) * (i+1));
		s->y = (float*)realloc(s->y, sizeof(float) * (i+1));
		s->d = (float*)realloc(s->d, sizeof(float) * (i+1));
		s->m = (float*)realloc(s->m, sizeof(float) * (i+1));
		fscanf(ent, "%f %f %f %f\n", &s->x[i], &s->y[i], &s->d[i], &s->m[i]);
		i++;
	}
	s->n = i;
	fclose(ent);
	system("rm -f _TMP");

//	fprintf(stderr, "END load %d.", s->n);
	return;
}

void reset(SET *s) {
	s->n = 0;
	s->x = NULL;
	s->y = NULL;
	s->m = NULL;
	s->d = NULL;

	s->y1 = 1992;
	s->y2 = 1992;
	s->m1 = 3;
	s->m2 = 9;
	s->d1 = 0;
	s->d2 = 1000;
	s->lat1 = -90;
	s->lat2 = 90;
	s->lon1 = -180;
	s->lon2 = 180;

	s->region = 0;

	return;
}

int findpt(SET *p, float ax, float ay) {

	int i;
	int current = -1;
	float d, temp;
	for(i = 0; i < p->n; i++) {
		temp = sqrt((ax - p->x[i])*(ax-p->x[i]) + (ay-p->y[i])* (ay-p->y[i]));
		if ( (current == -1) || (temp < d)) {
			d =temp;
			current = i;
		}
	}

	return current;
}

/*
 * Potting elements main functions
 */

int gomag(float *x, int n, float w, float **bins, float **freq) {
	float x1, x2;
	float xt;
	int i,j;

	/* Find x1, x2  as w multiple */
	minmax(x, n, &x1, &x2);

	for (xt = 0.0;xt <= x1; xt += w);
	x1 = xt - w;
	for (xt = 0.0 ;xt <= x2; xt += w);
	x2 = xt;

	/* Create bins */
	int nb = ((x2 - x1) / w);

	*bins = malloc(sizeof(float) * nb);
	*freq = malloc(sizeof(float) * nb);

	for(i=0; i<nb;  i++) (*bins)[i] = x1 + i*w;

	for(i=0; i < nb; i++) {
		float w1 = (*bins)[i] - w/2;
		(*freq)[i] = 0;
		for(j=0;j<n;j++)
			if (x[j] >= w1) (*freq)[i] += 1;
	}

	for(i=0; i < nb; i++) (*freq)[i] = log10((*freq)[i]);

	return nb;
}

int godep(float *x, int n, float w, float **bins, float **freq) {
	float x1, x2;
	float xt;
	int i,j;

	/* Find x1, x2  as w multiple */
	minmax(x, n, &x1, &x2);

	for (xt = 0.0;xt <= x1; xt += w);
	x1 = xt - w;
	for (xt = 0.0 ;xt <= x2; xt += w);
	x2 = xt;

	/* Create bins */
	int nb = ((x2 - x1) / w);

	*bins = malloc(sizeof(float) * nb);
	*freq = malloc(sizeof(float) * nb);

	for(i=0; i<nb;  i++) (*bins)[i] = x1 + i*w;

	for(i=0; i < nb; i++) {
		float w1 = (*bins)[i];
		float w2 = (*bins)[i] + w;
		(*freq)[i] = 0;
		for(j=0;j<n;j++)
			if (x[j] >= w1 && x[j] <= w2) (*freq)[i] += 1;
	}

	for(i=0; i < nb; i++) {
		if ((*freq)[i] > 0)
			(*freq)[i] = log10((*freq)[i]);
	}

	return nb;
}

void plothistogram(SET *p, float w, int mode) {
	float *x = (mode == MAG) ? p->m : p->d;
	int i;
	char t[1024];

	float x1, x2, y1, y2;

	float *bins = NULL;
	float *freq = NULL;
	float nb;

	float a, b;
	float rms = -1;

	if (mode == MAG) {
		nb = gomag(x, p->n, w, &bins, &freq);
		rms = linefit(bins, freq, nb, &a, &b);
	} else {
		nb = godep(x, p->n, w, &bins, &freq);
	}

	/*
	 * Plot
	 */
	cpgsvp(0.63, 0.93, 0.07, 0.30);

	minmax(x, p->n, &x1, &x2);
	minmax(freq, nb, &y1, &y2);
	cpgswin(x1, x2, y1, y2 * 1.2);
	cpgbox("BCNST", 0.0, 0, "BCMST", 0.0, 0);

	/*
	 * Labels
	 */
	cpgsch(0.7);
	cpgmtxt("L", 2.2, 0.0, 0.0, "[H] Trocar Mag/Dep");
	cpgmtxt("L", 1.0, 0.0, 0.0, "[B] Ajustar largura do bin");

	if (mode == MAG) {
		cpgmtxt("R", 3.0, 0.5, 0.5, "Log(n) Acumulado");
		cpgmtxt("B", 3.0, 0.5, 0.5, "Magnitude");
	} else {
		cpgmtxt("R", 3.0, 0.5, 0.5, "Log(n)");
		cpgmtxt("B", 3.0, 0.5, 0.5, "Profundidade (km)");
	}
	cpgsch(FS);

	/*
	 * Plots
	 */
	cpgbin(nb, bins, freq, 1);
	cpgmove(x1, a*x1 + b);
	cpgdraw(x2, a*x2 + b);

	if (mode == MAG) {
		cpgsch(0.7);
		sprintf(t,"f(x)=%.2f\\.x+%.2f",a,b);
		cpgmtxt("T",-2.0, 0.9, 1.0,t);

		sprintf(t,"b=%.2f",fabs(a));
		cpgmtxt("T",-3.2, 0.9, 1.0,t);
		cpgsch(FS);
	}

	cpgbbuf();

	/*
	 * Terminate
	 */
	cpgsci(1);
	cpgslw(1);

	cpgebuf();

	if (bins != NULL) free(bins);
	if (freq != NULL) free(freq);

	bins = NULL;
	freq = NULL;

	return;
}

void plotsection(SET *p, int mode) {
	float *x = (mode == LAT) ? p->y : p->x;
	float *y = p->d;
	int i;

	float x1, x2, y1, y2;

	x1 = x2 = x[0];
	y1 = y2 = y[0];
	for(i=0;i<p->n; i++) {
		if (x[i] < x1) x1 = x[i];
		if (x[i] > x2) x2 = x[i];
		if (y[i] < y1) y1 = y[i];
		if (y[i] > y2) y2 = y[i];
	}

	y2 = (( (int)y2 / 50 ) + 1) * 50.0;

	// Plot
	cpgsvp(0.07, 0.52, 0.07, 0.30);
	cpgswin(x1, x2, y2, 0.0);
	cpgbox("BCNST", 0.0, 0, "BCNST", 0.0, 0);
	cpgsch(0.7);
	cpgmtxt("R", 1.0, 0.0, 0.0, "[L] Trocar Lat/Lon");
	cpgmtxt("B", 3.0, 0.5, 0.5, (mode == LAT) ? "Latitude\\m94" : "Longitude\\m94");
	cpgmtxt("L", 3.0, 0.5, 0.5, "Profundidade (km)");
	cpgsch(FS);
	cpgbbuf();

	for(i = 0; i< p->n; i++) {
		cpgsci(depthcolor(y[i]));
		cpgpt1(x[i], y[i], 1);
	}
	cpgsci(1);


	// Terminate
	cpgsci(1);
	cpgslw(1);

	cpgebuf();
	return;
}

void plot(GRAPHCONTROL *gr, SET *p) {
	char t[1024];

	cpgsch(FS);
	cpgsvp(0.07, 0.93, 0.35, 0.9);
	cpgeras();
	cpgswin(gr->xmin, gr->xmax, gr->ymin, gr->ymax);
	cpgbox("BCNST", 0.0, 0, "BCNST", 0.0, 0);

	cpgbbuf();

	cpgsch(0.8);
	float yp = 3.4;

	sprintf(t,"[n] Ano: %d/%d", p->y1, p->y2);
	cpgmtxt("T", yp, 0.0, 0.0, t);

	sprintf(t,"[m] Magnitude: %.2f/%.2f", p->m1, p->m2);
	cpgmtxt("T", yp, 0.25, 0.0, t);

	sprintf(t,"[s/0] Selecionar Regiao");
	(p->region) ? cpgsci(ON) : cpgsci(OFF);
	cpgmtxt("T", yp, 0.6, 0.0, t);
	cpgsci(1);


	yp -= 1.2;

	sprintf(t,"N: %ld",p->n);
	cpgmtxt("T", yp, 0.0, 0.0, t);

	sprintf(t,"[p] Profundidade(p): %.1f/%.1f",p->d1, p->d2);
	cpgmtxt("T", yp, 0.25, 0.0, t);

	sprintf(t,"Longitude: %.2f/%.2f",p->lon1, p->lon2);
	cpgmtxt("T", yp, 0.6, 0.0, t);

	yp -= 1.2;

	sprintf(t,"Latitude: %.2f/%.2f", p->lat1, p->lat2);
	cpgmtxt("T", yp, 0.6, 0.0, t);

	sprintf(t,"[w] Zoom para todo o mapa");
	cpgmtxt("T", yp, 0.25, 0.0, t);



	sprintf(t,"[c] Cor: %s", (gr->colormode == COLORDEPTH) ? "Profundidade" : (gr->colormode == COLORMAG) ? "Magnitude" : "Neutra");
	cpgmtxt("R", 1.0, 1.0, 1.0, t);


	(gr->hascontinents) ? cpgsci(ON) : cpgsci(OFF);
	sprintf(t,"[1] Continentes");
	cpgmtxt("R", 1.0, 0.25, 0.0, t);

	(gr->hasplates) ? cpgsci(ON) : cpgsci(OFF);
	sprintf(t,"[2] Placas");
	cpgmtxt("R", 1.0, 0.0, 0.0, t);

	// Legenda cores
	cpgsci(1);
	cpgsch(FS);

		/* Graphs */
	int i;
	if (gr->haspoints) {
		if (gr->colormode == COLORDEPTH)
			for(i = 0; i< p->n; i++) {
				cpgsci(depthcolor(p->d[i]));
				cpgpt1(p->x[i], p->y[i], 1);
			}
		else if (gr->colormode == COLORMAG)
			for(i = 0; i< p->n; i++) {
				cpgsci(magcolor(p->m[i]));
				cpgpt1(p->x[i], p->y[i], 1);
			}
		else
			cpgpt(p->n, p->x, p->y, 1);
		cpgsci(1);
	}

	if (gr->hascontinents == 1) {
		cpgsci(1);
		cpgslw(2);

		for(i=0; i < ncontinentes; i++) {
			if (continentes[i][0] == -999 && continentes[i][1] == 999 ) {
				i++;
				cpgmove(continentes[i][0], continentes[i][1]);
				continue;
			}
			cpgdraw(continentes[i][0], continentes[i][1]);
		}
	}

	if (gr->hasplates == 1) {
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

	cpgebuf();
	return;
}

int enhance(GRAPHCONTROL *c, SET *p, int i) {
	float x, y, d, m;
	char t[1024];

	x = p->x[i];
	y = p->y[i];
	d = p->d[i];
	m = p->m[i];

	float xhalfwindow = (c->xmax + c->xmin) / 2;
	float yhalfwindow = (c->ymax + c->ymin) / 2;

			// Mark the dot
	cpgsch(FS);
	cpgsci(5);
	cpgpt1(x, y,-4);
	cpgsci(1);
	cpgsch(FS);

	// Info box
	if (x < xhalfwindow && y < yhalfwindow) {
		cpgsvp(0.55, 0.9, 0.7, 0.8);
	} else if (x < xhalfwindow && y > yhalfwindow) {
		cpgsvp(0.55, 0.9, 0.45, 0.55);
	} else if (x > xhalfwindow && y < yhalfwindow) {
		cpgsvp(0.12, 0.47, 0.7, 0.8);
	} else if (x > xhalfwindow && y > yhalfwindow) {
		cpgsvp(0.12, 0.47, 0.45, 0.55);
	}

	cpgswin(0.0, 1.0, 0.0, 1.0);
	cpgsci(1);
	cpgrect(0.0, 1.0, 0.0, 1.0);
	cpgsci(0);
	cpgrect(0.02, .98, 0.02, 0.98);
	cpgsci(1);

	sprintf(t,"Evento, %d",i);
	cpgmtxt("T", -1.0, 0.1, .0, t);
	sprintf(t,"Long. %.2f Lat. %.2f",x,y);
	cpgmtxt("T", -2.0, 0.1, .0, t);
	sprintf(t,"Prof. %.1f Mag. %.1f",d, m);
	cpgmtxt("T", -3.0, 0.1, .0, t);

	return -1;
}

/*
 * Main method
 */

int main(int argc, char **argv) {
	int ID;
	float x1, y1, x2, y2;
	char ch;

	int enhanceid = -1;
	int sectionmode = LON;
	int histmode = MAG;
	float binw = 0.25;

	reset(&mainset);
	load(&mainset);

	rangeadjust(&control, &mainset, NULL, NULL, NULL, NULL);

	ID = cpgopen("/xwindow");
	resizemax(0.8);
	cpgask(0);

	while (ch != 'Q') {
		plot(&control, &mainset);

		if (enhanceid != -1) {
			enhanceid = enhance(&control, &mainset, enhanceid);
		}

		if (mainset.region) {
			plotsection(&mainset, sectionmode);
			plothistogram(&mainset, binw, histmode);
		}

		// Restore default map position
		cpgsvp(0.07, 0.93, 0.35, 0.9);
		cpgswin(control.xmin, control.xmax, control.ymin, control.ymax);
		ch = getonechar(&x1, &y1, 0);
		switch(ch) {

		case('B'): {
			binw = lerfloat("Entre com o valor da largura do bin?");
			break;
		}

		case('A'): {
			enhanceid = findpt(&mainset, x1, y1);
			break;
		}

		case('2'): {
			control.hasplates = (control.hasplates == 1) ? 0 : 1;
			break;
		}

		case('1'): {
			control.hascontinents = (control.hascontinents == 1) ? 0 : 1;
			break;
		}

		case ('W'): {
			control.xmax = 180;
			control.xmin = -180;
			control.ymax = 90;
			control.ymin = -90;
			break;
		}

		case('0'): {
			mainset.lon1 = -180;
			mainset.lon2 = 180;
			mainset.lat1 = -90;
			mainset.lat2 = 90;
			load(&mainset);
			mainset.region = 0;
			break;
		}

		case ('L'): {
			sectionmode = (sectionmode == LAT) ? LON : LAT;
			break;
		}

		case ('H'): {
			histmode = (histmode == MAG) ? DEP : MAG;
			binw = (histmode == MAG) ? 0.25 : 25.0;
			break;
		}

		case ('S'): {
			x2 = x1;
			y2 = y1;
			getonechar(&x2, &y2, 1);
			order(&x1,&x2);
			order(&y1,&y2);
			mainset.lon1 = x1;
			mainset.lon2 = x2;
			mainset.lat1 = y1;
			mainset.lat2 = y2;
			load(&mainset);
			mainset.region = 1;
			break;
		}

		case('X'): { // ZOOM
			x2 = x1;
			y2 = y1;
			getonechar(&x2, &y2, 1);
			order(&x1,&x2);
			order(&y1,&y2);
			rangeadjust(&control, &mainset, &x1, &x2, &y1, &y2);
			break;
		}

		case('C'): { // Color mode
			control.colormode ++;
			if (control.colormode>COLORMAG) control.colormode = COLORNONE;
			break;
		}

		case('N'): { /* Ano */
			int y1, y2;
			y1 = lerint("Entre com ano inicial?");
			y2 = lerint("Entre com ano final?");
			orderint(&y1,&y2);
			mainset.y1 = y1;
			mainset.y2 = y2;
			load(&mainset);
			break;
		}

		case('M'): { /* Magnitude */
			float m1, m2;
			m1 = lerfloat("Entre com a Magnitude Inicial?");
			m2 = lerfloat("Entre com a Magnitude Final?");
			order(&m1,&m2);
			mainset.m1 = m1;
			mainset.m2 = m2;
			load(&mainset);
			break;
		}

		case('P'): { /* Profundidade */
			float d1, d2;
			d1 = lerfloat("Entre com a Profundidade Minima?");
			d2 = lerfloat("Entre com a Profundidade Maxima?");
			order(&d1,&d2);
			mainset.d1 = d1;
			mainset.d2 = d2;
			load(&mainset);
			break;
		}

		}
	}

	cpgclos();
}
