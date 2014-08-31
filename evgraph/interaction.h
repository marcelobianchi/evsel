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

#include <cpgplot.h>

#define ERROR 2
#define WARN 1
#define OK 0

/* Interaction */
extern char message[2048];

void alert(int level);
float lerfloat(char *text);
void lerchar(char *text, char *output, int max);
int lerint(char *text);
int getmouse(int a, int b, int c, int d, char *message);
char getonechar(float *axx, float *ayy,  int pos);
int yesno(char *text);
