#include <stdio.h>
//#include <windows.h>	   // Standard header for MS Windows applications
#include <GL/gl.h>		   // Open Graphics Library (OpenGL) header
#include <GL/glut.h>
#include <math.h>

#ifndef CALLBACK
#define CALLBACK
#endif

#define PI 3.14
#define NONE -1
#define POINT_SIZE 10
#define CHECKED 1
#define UNCHECKED 0
#define SIZE 1.0
static GLfloat line_red=0.0;
static GLfloat line_green=0.0;
static GLfloat line_blue=0.0;

static GLfloat fill_red=1.0;
static GLfloat fill_green=1.0;
static GLfloat fill_blue=1.0;

GLdouble currentWinding = GLU_TESS_WINDING_NONZERO;
GLUtesselator *tobj;
GLuint list;

void CALLBACK beginCallback(GLenum which);
void CALLBACK endCallback();
void CALLBACK errorCallback(GLenum errorCode);
void paintPolygon(int i,int z);

typedef struct edge{
    int x1,x2,y1,y2;
    int cnt = 0;
    double txy[20][3];
}edge;

typedef struct polygon{
    edge polygonEdges[20]; //polygon edges
    int cnt=0; //how many polygon edges
    int w_num = -1000;      //winding number
    int isClockWise;
    int visited = 0;
    int isCounterClockWise;
}polygon;

typedef struct dividedEdge{
	int x1,x2,y1,y2;
	int normal=UNCHECKED;
	int backward=UNCHECKED;
}dividedEdge;


double alpha=600.0;
double beta=600.0;
double angle =0.0;
double angle2 = 1.0;

polygon polygonList[40];
static int external = -1;  //position of the external polygon in the polygon list.
int externalPolygons[30];
int externalPolygonsCnt = 0;
static int polygoncnt = 0;
static int option,button,x,y,color_flag=0;
static int xi,yi;
static int found = -1;
static int flag_normal=0;
static int flag_backward=0;
static int window;
static int cordpnt=-1;
static int edgecnt = 0;
static int foundcnt = 0;
static int cord[20][2]; //20 shmeia x,y
edge edgeList[40];  //oles oi akmes
dividedEdge edgeFound[10];
dividedEdge dividedEdgeList[80];
int polygonsToExtrude[40];
int pol = 0;
int xCircle,zCircle;
int extrudeType = -1;
int rule =- 1;
int firstTimeOpened = 0;
void display(void);
void createMenu(void);
void addToList(double t1,double t2,int xpr,int ypr,int x1,int x2,int x3,int x4,int y1,int y2,int y3,int y4,int edge1,int edge2);
void inter(int x1, int y1, int x2, int y2, int x3, int y3, int x4, int y4, int edge1, int edge2);
void pickPolygonsToFill();


void Init();

void CALLBACK beginCallback(GLenum which)
{
   glBegin(which);
}

void CALLBACK errorCallback(GLenum errorCode)
{
   const GLubyte *estring;

   estring = gluErrorString(errorCode);
   fprintf(stderr, "Tessellation Error: %s\n", estring);
   exit(0);
}

void CALLBACK endCallback(void)
{
   glEnd();
}


void findIntersectionPoints2(){
	int i,j, x1, y1, x2, y2, x3, y3, x4, y4;
	for(i=0;i<(cordpnt+1);i++){
		for(j=i+1;j<(cordpnt+1);j++){
			//edge 1
			x1=cord[i][0];
			y1=cord[i][1];
			if(i==cordpnt){
				x2=cord[0][0];
				y2=cord[0][1];
			} else {
				x2=cord[i+1][0];
				y2=cord[i+1][1];
			}

			//edge 2
			x3=cord[j][0];
			y3=cord[j][1];
			if(j==cordpnt){
				x4=cord[0][0];
				y4=cord[0][1];
			} else {
				x4=cord[j+1][0];
				y4=cord[j+1][1];
			}

			if(i!=j){
				inter(x1, y1, x2, y2, x3, y3, x4, y4, i, j);
			}
		}
	}
}
void inter(int x1, int y1, int x2, int y2, int x3, int y3, int x4, int y4, int edge1, int edge2){
	int xpr, ypr;
	double t1,t2;
	// Line AB represented as a1x + b1y = c1
    double a1 = y2 - y1;
    double b1 = x1 - x2;
    double c1 = a1*(x1) + b1*(y1);

    // Line CD represented as a2x + b2y = c2
    double a2 = y4 - y3;
    double b2 = x3 - x4;
    double c2 = a2*(x3)+ b2*(y3);

    double determinant = a1*b2 - a2*b1;

    if (determinant == 0)
    {
        // The lines are parallel. This is simplified
        // by returning a pair of FLT_MAX
    }
    else
    {
        double x = (b2*c1 - b1*c2)/determinant;
        double y = (a1*c2 - a2*c1)/determinant;
        xpr=(int)(x-0.5);
	ypr=(int)(y-0.5);

		//AFOU VREI MIA TOMH, VRISKEI TO t THS STIS DUO EUTHEIES KAI KALEI ADDTOLIST
		if((x2-x1)!=0)
			t1 =((x - x1)/double(x2-x1));
		else{
			t1 = (y-y1)/double(y2-y1);
		}
		if((x3-x4)!=0){
			t2 = (x - x3)/(double)(x4-x3);
		}else{
			t2 = (y - y3)/(double)(y4-y3);
		}

		if((t1>0 && t1<1) && (t2>0 && t2<1)){
			addToList(t1,t2,xpr,ypr,x1,x2,x3,x4,y1,y2,y3,y4,edge1,edge2);
		}



    }

}


void addToList(double t1,double t2,int xpr,int ypr,int x1,int x2,int x3,int x4,int y1,int y2,int y3,int y4,int edge1,int edge2){
	edgeList[edge1].x1=x1;
	edgeList[edge1].y1=y1;
	edgeList[edge1].x2=x2;
	edgeList[edge1].y2=y2;
	edgeList[edge1].txy[edgeList[edge1].cnt][0]=t1;
	edgeList[edge1].txy[edgeList[edge1].cnt][1]=xpr;
	edgeList[edge1].txy[edgeList[edge1].cnt][2]=ypr;
	edgeList[edge1].cnt++;

	edgeList[edge2].x1=x3;
	edgeList[edge2].y1=y3;
	edgeList[edge2].x2=x4;
	edgeList[edge2].y2=y4;
	edgeList[edge2].txy[edgeList[edge2].cnt][0]=t2;
	edgeList[edge2].txy[edgeList[edge2].cnt][1]=xpr;
	edgeList[edge2].txy[edgeList[edge2].cnt][2]=ypr;
	edgeList[edge2].cnt++;


}

void swap(double *xp, double *yp)
{
    double temp = *xp;
    *xp = *yp;
    *yp = temp;
}

void selectionSort(double arr[40][3], int n)
{
    int i, j, min_idx;

    // One by one move boundary of unsorted subarray
    for (i = 0; i < n-1; i++)
    {
        // Find the minimum element in unsorted array
        min_idx = i;
        for (j = i+1; j < n; j++)
          if (arr[j][0] < arr[min_idx][0])
            min_idx = j;

        // Swap the found minimum element with the first element
        swap(&arr[min_idx][0], &arr[i][0]);
        swap(&arr[min_idx][1], &arr[i][1]);
        swap(&arr[min_idx][2], &arr[i][2]);
    }
}


void divideEdges(){
	int i,j;
	edgecnt = 0;
	for(i=0;i<(cordpnt+1);i++){
		selectionSort(edgeList[i].txy,edgeList[i].cnt);
		
	}
	for(i=0;i<(cordpnt+1);i++){
		for (j=0;j<edgeList[i].cnt;j++){
			if(j==0){
				dividedEdgeList[edgecnt].x1 = edgeList[i].x1;
				dividedEdgeList[edgecnt].y1 = edgeList[i].y1;
				dividedEdgeList[edgecnt].x2 = edgeList[i].txy[j][1];
				dividedEdgeList[edgecnt].y2 = edgeList[i].txy[j][2];
				dividedEdgeList[edgecnt].normal=UNCHECKED;
	 			dividedEdgeList[edgecnt].backward=UNCHECKED;
			}
			else{
				dividedEdgeList[edgecnt].x1 = edgeList[i].txy[j-1][1];
				dividedEdgeList[edgecnt].y1 = edgeList[i].txy[j-1][2];
				dividedEdgeList[edgecnt].x2 = edgeList[i].txy[j][1];
				dividedEdgeList[edgecnt].y2 = edgeList[i].txy[j][2];
				dividedEdgeList[edgecnt].normal=UNCHECKED;
	 			dividedEdgeList[edgecnt].backward=UNCHECKED;
			}
			edgecnt++;
		}
		if(edgeList[i].cnt > 0){
			dividedEdgeList[edgecnt].x1 = edgeList[i].txy[edgeList[i].cnt-1][1];
			dividedEdgeList[edgecnt].y1 = edgeList[i].txy[edgeList[i].cnt-1][2];
			dividedEdgeList[edgecnt].x2 = edgeList[i].x2;
			dividedEdgeList[edgecnt].y2 = edgeList[i].y2;
			dividedEdgeList[edgecnt].normal=UNCHECKED;
	 		dividedEdgeList[edgecnt].backward=UNCHECKED;

		}
		else{
			dividedEdgeList[edgecnt].x1 = edgeList[i].x1;
			dividedEdgeList[edgecnt].y1 = edgeList[i].y1;
			dividedEdgeList[edgecnt].x2 = edgeList[i].x2;
			dividedEdgeList[edgecnt].y2 = edgeList[i].y2;
			dividedEdgeList[edgecnt].normal=UNCHECKED;
	 		dividedEdgeList[edgecnt].backward=UNCHECKED;
		}
		edgecnt++;
	}

	
}


int findEdge(int x_start,int y_start,int x,int y){
	int i,mini,maxi;
	double min,max;
	int det[20],detcnt = 0,ab,flag = 0;
	double angle[20],dist_a,dist_b,c;
	flag_normal=0;
	flag_backward=0;
	foundcnt = 0;
	for(i=0;i<edgecnt;i++) {
		if(dividedEdgeList[i].x1==x && dividedEdgeList[i].y1==y){
			if(dividedEdgeList[i].x2==x_start && dividedEdgeList[i].y2==y_start){
				edgeFound[foundcnt].x1=x_start;
				edgeFound[foundcnt].y1=y_start;
				edgeFound[foundcnt].x2=x;
				edgeFound[foundcnt].y2=y;
			}else{
				edgeFound[foundcnt].x1=x;
				edgeFound[foundcnt].y1=y;
				edgeFound[foundcnt].x2=dividedEdgeList[i].x2;
				edgeFound[foundcnt].y2=dividedEdgeList[i].y2;
			}
			foundcnt++;
		}else if(dividedEdgeList[i].x2==x && dividedEdgeList[i].y2==y){
			if(dividedEdgeList[i].x1==x_start && dividedEdgeList[i].y1==y_start){
				edgeFound[foundcnt].x1=x_start;
				edgeFound[foundcnt].y1=y_start;
				edgeFound[foundcnt].x2=x;
				edgeFound[foundcnt].y2=y;
			}else{
				edgeFound[foundcnt].x1=x;
				edgeFound[foundcnt].y1=y;
				edgeFound[foundcnt].x2=dividedEdgeList[i].x1;
				edgeFound[foundcnt].y2=dividedEdgeList[i].y1;
			}
			foundcnt++;

		}
	}
	for(i=0;i<foundcnt;i++){
		det[i] = (edgeFound[i].x2 - edgeFound[i].x1) * (y - y_start) - (edgeFound[i].y2 - edgeFound[i].y1) * (x - x_start);
		ab = (edgeFound[i].x2 - edgeFound[i].x1) * (x - x_start) + (edgeFound[i].y2 - edgeFound[i].y1)* (y - y_start);
		dist_a=sqrt((double)(pow((double)(edgeFound[i].x2 - edgeFound[i].x1),2.0) + pow((double)(edgeFound[i].y2 - edgeFound[i].y1),2.0)));
		dist_b= sqrt((double)(pow((double)(x - x_start),2.0) + pow((double)(y-y_start),2.0)));
		c = ab/(double)(dist_a*dist_b);
		angle[i] = acos(c);
	}

	min = 3.15;
	max = 0;
	mini = 0;
	maxi = 0;
	for(i=0;i<foundcnt;i++){
		if(det[i] > 0){
			flag = 1;
		}
		if(flag == 1 && det[i] > 0 && angle[i] > max){
			max = angle[i];
			maxi = i;
		}
		else if(flag == 0 && det[i] < 0 && angle[i] < min){
			min = angle[i];
			mini = i;
		}
	}

	if(flag == 1){
		for(i=0;i<edgecnt;i++){
			if((dividedEdgeList[i].x1 == edgeFound[maxi].x1 && dividedEdgeList[i].x2 == edgeFound[maxi].x2 ) &&
				(dividedEdgeList[i].y1 == edgeFound[maxi].y1 && dividedEdgeList[i].y2 == edgeFound[maxi].y2)){
				flag_normal=1;
				return i;
			}else if(dividedEdgeList[i].x1 == edgeFound[maxi].x2 && dividedEdgeList[i].x2 == edgeFound[maxi].x1 &&
				dividedEdgeList[i].y1 == edgeFound[maxi].y2 && dividedEdgeList[i].y2 == edgeFound[maxi].y1){
				flag_backward=1;
				return i;
			}
		}
	}
	else{
		for(i=0;i<edgecnt;i++){
			if((dividedEdgeList[i].x1 == edgeFound[mini].x1 && dividedEdgeList[i].x2 == edgeFound[mini].x2) && 			//normal
				(dividedEdgeList[i].y1 == edgeFound[mini].y1 && dividedEdgeList[i].y2 == edgeFound[mini].y2)){
				//dividedEdge[i].normal=CHECKED;
				flag_normal=1;
				return i;
			}else if((dividedEdgeList[i].x1 == edgeFound[mini].x2 && dividedEdgeList[i].x2 == edgeFound[mini].x1 ) &&	//backwards
				(dividedEdgeList[i].y1 == edgeFound[mini].y2 && dividedEdgeList[i].y2 == edgeFound[mini].y1)){
				//dividedEdge[i].backward=CHECKED;
				flag_backward=1;
				return i;
			}
		}

	}
	return -1;
}

void findPolygons(){
	int i,j,xs,ys,xe,ye,loop=1,k,sum = 0;
	foundcnt=0;
	for(i=0;i<polygoncnt;i++){
        polygonList[i].cnt = 0;
	}
	polygoncnt = 0;
	for(i=0;i<edgecnt;i++){
		if(dividedEdgeList[i].normal==UNCHECKED){
			xs=dividedEdgeList[i].x1;
			ys=dividedEdgeList[i].y1;
			xe=dividedEdgeList[i].x2;
			ye=dividedEdgeList[i].y2;
			dividedEdgeList[i].normal = CHECKED;
			polygonList[polygoncnt].polygonEdges[polygonList[polygoncnt].cnt].x1 = xs;  //add first polygon edge to the first polygon.
			polygonList[polygoncnt].polygonEdges[polygonList[polygoncnt].cnt].y1 = ys;
			polygonList[polygoncnt].polygonEdges[polygonList[polygoncnt].cnt].x2 = xe;
			polygonList[polygoncnt].polygonEdges[polygonList[polygoncnt].cnt].y2 = ye;
			polygonList[polygoncnt].cnt ++;
			while(1){
				k = findEdge(xs,ys,xe,ye);
				if(dividedEdgeList[k].normal==UNCHECKED && flag_normal==1){
					dividedEdgeList[k].normal=CHECKED;
					xs=dividedEdgeList[k].x1;
					ys=dividedEdgeList[k].y1;
					xe=dividedEdgeList[k].x2;
					ye=dividedEdgeList[k].y2;
                    polygonList[polygoncnt].polygonEdges[polygonList[polygoncnt].cnt].x1 = xs;  //add all polygon edges to the polygon.
                    polygonList[polygoncnt].polygonEdges[polygonList[polygoncnt].cnt].y1 = ys;
                    polygonList[polygoncnt].polygonEdges[polygonList[polygoncnt].cnt].x2 = xe;
                    polygonList[polygoncnt].polygonEdges[polygonList[polygoncnt].cnt].y2 = ye;
                    polygonList[polygoncnt].cnt ++;
				}
				if(dividedEdgeList[k].backward==UNCHECKED && flag_backward==1){
					dividedEdgeList[k].backward=CHECKED;
					xs=dividedEdgeList[k].x2;
					ys=dividedEdgeList[k].y2;
					xe=dividedEdgeList[k].x1;
					ye=dividedEdgeList[k].y1;
                    polygonList[polygoncnt].polygonEdges[polygonList[polygoncnt].cnt].x1 = xs;  //add all polygon edges to the polygon.
                    polygonList[polygoncnt].polygonEdges[polygonList[polygoncnt].cnt].y1 = ys;
                    polygonList[polygoncnt].polygonEdges[polygonList[polygoncnt].cnt].x2 = xe;
                    polygonList[polygoncnt].polygonEdges[polygonList[polygoncnt].cnt].y2 = ye;
                    polygonList[polygoncnt].cnt ++;

				}
				if(dividedEdgeList[k].normal==CHECKED && dividedEdgeList[k].x1==dividedEdgeList[i].x1 && dividedEdgeList[k].y1==dividedEdgeList[i].y1
					&& dividedEdgeList[k].x2==dividedEdgeList[i].x2 && dividedEdgeList[k].y2==dividedEdgeList[i].y2){
                    polygonList[polygoncnt].w_num=-1000;
		    		polygonList[polygoncnt].visited = 0;
                    polygoncnt ++;  //point to the next thesis to add a new polygon
					break; //found polygon
				}

			}
		}

		if(dividedEdgeList[i].backward==UNCHECKED){
			xs=dividedEdgeList[i].x2;
			ys=dividedEdgeList[i].y2;
			xe=dividedEdgeList[i].x1;
			ye=dividedEdgeList[i].y1;
			dividedEdgeList[i].backward = CHECKED;
			polygonList[polygoncnt].polygonEdges[polygonList[polygoncnt].cnt].x1 = xs;  //add first polygon edge to the first polygon.
			polygonList[polygoncnt].polygonEdges[polygonList[polygoncnt].cnt].y1 = ys;
			polygonList[polygoncnt].polygonEdges[polygonList[polygoncnt].cnt].x2 = xe;
			polygonList[polygoncnt].polygonEdges[polygonList[polygoncnt].cnt].y2 = ye;
			polygonList[polygoncnt].cnt ++;
			while(1){
				k = findEdge(xs,ys,xe,ye);
				if(dividedEdgeList[k].normal==UNCHECKED && flag_normal==1){
					dividedEdgeList[k].normal=CHECKED;
					xs=dividedEdgeList[k].x1;
					ys=dividedEdgeList[k].y1;
					xe=dividedEdgeList[k].x2;
					ye=dividedEdgeList[k].y2;
					polygonList[polygoncnt].polygonEdges[polygonList[polygoncnt].cnt].x1 = xs;  //add all polygon edges to the polygon.
                    polygonList[polygoncnt].polygonEdges[polygonList[polygoncnt].cnt].y1 = ys;
                    polygonList[polygoncnt].polygonEdges[polygonList[polygoncnt].cnt].x2 = xe;
                    polygonList[polygoncnt].polygonEdges[polygonList[polygoncnt].cnt].y2 = ye;
                    polygonList[polygoncnt].cnt ++;
				}
				if(dividedEdgeList[k].backward==UNCHECKED && flag_backward==1){
					dividedEdgeList[k].backward=CHECKED;
					xs=dividedEdgeList[k].x2;
					ys=dividedEdgeList[k].y2;
					xe=dividedEdgeList[k].x1;
					ye=dividedEdgeList[k].y1;
					polygonList[polygoncnt].polygonEdges[polygonList[polygoncnt].cnt].x1 = xs;  //add all polygon edges to the polygon.
                    polygonList[polygoncnt].polygonEdges[polygonList[polygoncnt].cnt].y1 = ys;
                    polygonList[polygoncnt].polygonEdges[polygonList[polygoncnt].cnt].x2 = xe;
                    polygonList[polygoncnt].polygonEdges[polygonList[polygoncnt].cnt].y2 = ye;
                    polygonList[polygoncnt].cnt ++;
				}
				if(dividedEdgeList[k].backward==CHECKED && dividedEdgeList[k].x1==dividedEdgeList[i].x1 && dividedEdgeList[k].y1==dividedEdgeList[i].y1
					&& dividedEdgeList[k].x2==dividedEdgeList[i].x2 && dividedEdgeList[k].y2==dividedEdgeList[i].y2){
                    polygonList[polygoncnt].w_num=-1000;
		    polygonList[polygoncnt].visited = 0;
                    polygoncnt ++;  //point to the next thesis to add a new polygon
					break; //found polygon
				}
			}
		}
	}
}

int mostLeftEdge(){
	int i,j,minEcnt,twoMinI[2],out=-1;
	double ab,distA,angle[2],det[2];
	dividedEdge min,twoMin[2];
	int flag=0;
	min.x1=600;
	min.y1=600;
	min.x2=500;
	min.y2=500;
	double maxAngle=0.0;
	double maxAngle2=0.0;
	for(i=0;i<edgecnt;i++){
		if(dividedEdgeList[i].y1!=dividedEdgeList[i].y2){
			if((dividedEdgeList[i].x1 +  dividedEdgeList[i].y1 ) < (min.x1 + min.y1)||
				(dividedEdgeList[i].x2 + dividedEdgeList[i].y2 )< (min.x1 + min.y1)){
				min.x1=dividedEdgeList[i].x1;
				min.y1=dividedEdgeList[i].y1;
				min.x2=dividedEdgeList[i].x2;
				min.y2=dividedEdgeList[i].y2;
			}
		}
	}
	minEcnt=0;
	
	for(i=0;i<edgecnt;i++){
		if(dividedEdgeList[i].x1==min.x1 && dividedEdgeList[i].y1==min.y1){
			twoMin[minEcnt].x1=dividedEdgeList[i].x1;
			twoMin[minEcnt].y1=dividedEdgeList[i].y1;
			twoMin[minEcnt].x2=dividedEdgeList[i].x2;
			twoMin[minEcnt].y2=dividedEdgeList[i].y2;
			twoMinI[minEcnt]=i;
			minEcnt++;

		}else if(dividedEdgeList[i].x2==min.x1 && dividedEdgeList[i].y2==min.y1){
			twoMin[minEcnt].x1=dividedEdgeList[i].x2;
			twoMin[minEcnt].y1=dividedEdgeList[i].y2;
			twoMin[minEcnt].x2=dividedEdgeList[i].x1;
			twoMin[minEcnt].y2=dividedEdgeList[i].y1;		
			twoMinI[minEcnt]=i;
			minEcnt++;
		}
	}
	for(j=0;j<minEcnt;j++){
		ab=(twoMin[j].x2-twoMin[j].x1);
		distA=sqrt((double)(pow((double)(twoMin[j].x2-twoMin[j].x1),2.0) + pow((double)(twoMin[j].y2-twoMin[j].y1),2.0)));
		angle[j]=acos(ab/(double)distA);
		det[j]=-(twoMin[j].y2-twoMin[j].y1);
		
	}



	maxAngle=0.0;
	maxAngle2=0.0;
	for(i=0;i<minEcnt;i++){
		if(det[i] < 0){
			flag = 1;
		}
		if(flag == 1 && det[i] < 0 && angle[i] > maxAngle){
			maxAngle=angle[i];
			out=twoMinI[i];
		}
		else if(flag == 0 && det[i] > 0 && angle[i] > maxAngle2){
			maxAngle2=angle[i];
			out=twoMinI[i];
		}
	}

	return out;
}

void findExternalPolygons(){
	int i,j,k,flag;
	externalPolygonsCnt=0;
	for(i=0;i<polygoncnt;i++){
		if(i != external){
			for(j=0;j<polygonList[i].cnt;j++){
				for(k=0;k<polygonList[external].cnt;k++){
					if(polygonList[i].polygonEdges[j].x1 == polygonList[external].polygonEdges[k].x2 &&
						polygonList[i].polygonEdges[j].y1 == polygonList[external].polygonEdges[k].y2 &&
						polygonList[i].polygonEdges[j].x2 == polygonList[external].polygonEdges[k].x1 &&
						polygonList[i].polygonEdges[j].y2 == polygonList[external].polygonEdges[k].y1 ){

							externalPolygons[externalPolygonsCnt] = i;
							externalPolygonsCnt++;
							flag = 1;
							break;
					}
				}
				if(flag == 1){
					flag = 0;
					break;
				}
			}
		}
	}
}



void findExternalWindings(){
	int i,j,z,k,minEcnt,twoMinI[2],out=-1;
	int flag=0,isExtEdge = 0;
	double ab,distA,angle[2],det[2];
	dividedEdge min,twoMin[2];
	double maxAngle,maxAngle2;
	for(i=0;i<externalPolygonsCnt;i++){
		min.x1=600;
		min.y1=600;
		min.x2=500;
		min.y2=500;
		maxAngle=0.0;
		isExtEdge=0;

		for(j=0;j<polygonList[externalPolygons[i]].cnt;j++){
			if(polygonList[externalPolygons[i]].polygonEdges[j].y1!=polygonList[externalPolygons[i]].polygonEdges[j].y2){
				if((polygonList[externalPolygons[i]].polygonEdges[j].x1 +  polygonList[externalPolygons[i]].polygonEdges[j].y1 ) < (min.x1 + min.y1)||
					(polygonList[externalPolygons[i]].polygonEdges[j].x2 + polygonList[externalPolygons[i]].polygonEdges[j].y2 )< (min.x1 + min.y1)){
					min.x1=polygonList[externalPolygons[i]].polygonEdges[j].x1;
					min.y1=polygonList[externalPolygons[i]].polygonEdges[j].y1;
					min.x2=polygonList[externalPolygons[i]].polygonEdges[j].x2;
					min.y2=polygonList[externalPolygons[i]].polygonEdges[j].y2;
				}
			}
		}
		minEcnt=0;
			
		for(z=0;z<polygonList[externalPolygons[i]].cnt;z++){
			if(polygonList[externalPolygons[i]].polygonEdges[z].x1==min.x1 && polygonList[externalPolygons[i]].polygonEdges[z].y1==min.y1){
				twoMin[minEcnt].x1=polygonList[externalPolygons[i]].polygonEdges[z].x1;
				twoMin[minEcnt].y1=polygonList[externalPolygons[i]].polygonEdges[z].y1;
				twoMin[minEcnt].x2=polygonList[externalPolygons[i]].polygonEdges[z].x2;
				twoMin[minEcnt].y2=polygonList[externalPolygons[i]].polygonEdges[z].y2;
				
				twoMinI[minEcnt]=z;
				minEcnt++;

			}else if(polygonList[externalPolygons[i]].polygonEdges[z].x2==min.x1 && polygonList[externalPolygons[i]].polygonEdges[z].y2==min.y1){
				twoMin[minEcnt].x1=polygonList[externalPolygons[i]].polygonEdges[z].x2;
				twoMin[minEcnt].y1=polygonList[externalPolygons[i]].polygonEdges[z].y2;
				twoMin[minEcnt].x2=polygonList[externalPolygons[i]].polygonEdges[z].x1;
				twoMin[minEcnt].y2=polygonList[externalPolygons[i]].polygonEdges[z].y1;
					
				twoMinI[minEcnt]=z;
				minEcnt++;
			}
		}
		for(k=0;k<minEcnt;k++){
			ab=(twoMin[k].x2-twoMin[k].x1);
			distA=sqrt((double)(pow((double)(twoMin[k].x2-twoMin[k].x1),2.0) + pow((double)(twoMin[k].y2-twoMin[k].y1),2.0)));
			angle[k]=acos(ab/(double)distA);
			det[k]=-(twoMin[k].y2-twoMin[k].y1);
		}

		maxAngle=0.0;
		maxAngle2=0.0;
		for(k=0;k<minEcnt;k++){
			if(det[k] < 0){
				flag = 1;
			}
			if(flag == 1 && det[k] < 0 && angle[k] > maxAngle){
				maxAngle=angle[k];
				out=twoMinI[k];
			}
			else if(flag == 0 && det[k] > 0 && angle[k] > maxAngle2){
				maxAngle2=angle[k];
				out=twoMinI[k];
			}
		}

		for(j=0;j<polygonList[external].cnt;j++){
			if(polygonList[externalPolygons[i]].polygonEdges[out].x1 == polygonList[external].polygonEdges[j].x2 &&
				polygonList[externalPolygons[i]].polygonEdges[out].y1 == polygonList[external].polygonEdges[j].y2 &&
				polygonList[externalPolygons[i]].polygonEdges[out].x2 == polygonList[external].polygonEdges[j].x1 &&
				polygonList[externalPolygons[i]].polygonEdges[out].y2 == polygonList[external].polygonEdges[j].y1){
				isExtEdge = 1;
			}
		}
		if(isExtEdge == 0 && out == twoMinI[0]){
			out = twoMinI[1];
		}
		else if(isExtEdge == 0 && out == twoMinI[1]){
			out = twoMinI[0];
		}

		for(j=0;j<edgecnt;j++){
			if(((polygonList[externalPolygons[i]].polygonEdges[out].x1 == dividedEdgeList[j].x1) && 
			(polygonList[externalPolygons[i]].polygonEdges[out].y1 == dividedEdgeList[j].y1)  &&
			(polygonList[externalPolygons[i]].polygonEdges[out].x2 == dividedEdgeList[j].x2)  &&
			(polygonList[externalPolygons[i]].polygonEdges[out].y2 == dividedEdgeList[j].y2)) ||
			((polygonList[externalPolygons[i]].polygonEdges[out].x2 == dividedEdgeList[j].x1) &&
			(polygonList[externalPolygons[i]].polygonEdges[out].y2 == dividedEdgeList[j].y1) &&
			(polygonList[externalPolygons[i]].polygonEdges[out].x1 == dividedEdgeList[j].x2) &&
			(polygonList[externalPolygons[i]].polygonEdges[out].y1 == dividedEdgeList[j].y2))) {

				if(isExtEdge==0){
					if(dividedEdgeList[j].x1 > dividedEdgeList[j].x2){   //going left winding =-1
						polygonList[externalPolygons[i]].w_num = -1;
					}else if(dividedEdgeList[j].x1 < dividedEdgeList[j].x2){ //going right winding =1
						polygonList[externalPolygons[i]].w_num = 1;
					}

				}else
					if(dividedEdgeList[j].y1 > dividedEdgeList[j].y2){   //going down winding =1
						polygonList[externalPolygons[i]].w_num = 1;
					}else if(dividedEdgeList[j].y1 < dividedEdgeList[j].y2){ //going up winding =-1
						polygonList[externalPolygons[i]].w_num = -1;
					}
			}
		}

	}
	
	return;
}


void findInternalWindingNumbers(){
	int i,j,h,k,m,n,b,flag = 0,internalcnt=1,polygons[40],count=-1,flag2=0,currentPoly;

	for(i=0;i<externalPolygonsCnt;i++){			//arxikopoihsh tou polygons sta ekswterika.. tha ksekinisei apta ekswterika
		count++; 								//kai tha proxwrhsei pros ta eswterika sthn while..
		polygons[count] = externalPolygons[i];
		polygonList[polygons[count]].visited = 1;
	}
	internalcnt=polygoncnt-count-1;
	currentPoly = polygonList[polygons[count]].cnt;	//starting from the end of the list

	while(count>=0){
		for(j=0;j<currentPoly;j++){
			for(k=0;k<polygoncnt;k++){
				if(k != external && k != polygons[count]){
					for(h=0;h<polygonList[k].cnt;h++){
						if(polygonList[polygons[count]].polygonEdges[j].x1 == polygonList[k].polygonEdges[h].x2 &&
							polygonList[polygons[count]].polygonEdges[j].y1 == polygonList[k].polygonEdges[h].y2 &&
							polygonList[polygons[count]].polygonEdges[j].x2 == polygonList[k].polygonEdges[h].x1 &&
							polygonList[polygons[count]].polygonEdges[j].y2 == polygonList[k].polygonEdges[h].y1 &&
							polygonList[k].visited == 0){

							for(b=0;b<edgecnt;b++){
								if(((polygonList[polygons[count]].polygonEdges[j].x1 == dividedEdgeList[b].x1) && 
								(polygonList[polygons[count]].polygonEdges[j].y1 == dividedEdgeList[b].y1)  &&
								(polygonList[polygons[count]].polygonEdges[j].x2 == dividedEdgeList[b].x2)  &&
								(polygonList[polygons[count]].polygonEdges[j].y2 == dividedEdgeList[b].y2))){
									polygonList[k].w_num = polygonList[polygons[count]].w_num + 1;
									polygons[count] = k;
									flag2 = 1;
									polygonList[k].visited = 1;
									break;
								}
								else if (((polygonList[polygons[count]].polygonEdges[j].x2 == dividedEdgeList[b].x1) && 
								(polygonList[polygons[count]].polygonEdges[j].y2 == dividedEdgeList[b].y1)  &&
								(polygonList[polygons[count]].polygonEdges[j].x1 == dividedEdgeList[b].x2)  &&
								(polygonList[polygons[count]].polygonEdges[j].y1 == dividedEdgeList[b].y2))){
									polygonList[k].w_num = polygonList[polygons[count]].w_num - 1;
									polygons[count] = k;
									flag2 = 1;
									polygonList[k].visited = 1;
									break;
								}
							}
							if(flag2 == 1){
								break;
							}
						}
					}
					if(flag2 == 1){
						break;
					}
				}
			}
		}
		if(flag2 == 1){
			flag2 = 0;
		}
		else{
			count--;
		}
		currentPoly = polygonList[polygons[count]].cnt;
	}
	
}


void findWindingNumbers(){
	int mostLeft,i,j;
	int goingDown=0;
	int goingUp=0;
	mostLeft=mostLeftEdge();
	if(dividedEdgeList[mostLeft].y1>dividedEdgeList[mostLeft].y2){
		goingDown=1;
	}else{
		goingUp=1;
	}
	
	for(i=0;i<polygoncnt;i++){
        for(j=0;j<polygonList[i].cnt;j++){
			if(goingUp==1){
				if(dividedEdgeList[mostLeft].x2 == polygonList[i].polygonEdges[j].x1 && dividedEdgeList[mostLeft].y2 == polygonList[i].polygonEdges[j].y1 &&
				dividedEdgeList[mostLeft].x1 == polygonList[i].polygonEdges[j].x2 && dividedEdgeList[mostLeft].y1 == polygonList[i].polygonEdges[j].y2 ) {
					external = i;//Find the i of the external polygon
				}
			}
			if(goingDown == 1){
				if(dividedEdgeList[mostLeft].x1 == polygonList[i].polygonEdges[j].x1 && dividedEdgeList[mostLeft].y1 == polygonList[i].polygonEdges[j].y1 &&
				dividedEdgeList[mostLeft].x2 == polygonList[i].polygonEdges[j].x2 && dividedEdgeList[mostLeft].y2 == polygonList[i].polygonEdges[j].y2 ) {
					external = i;
				}
       		}
		}
    }

    findExternalPolygons();
    findExternalWindings();
    findInternalWindingNumbers();
}

void pickPolygonsToFill(){
	int i;
        pol = 0;
	if(rule == 0){
		for(i=0;i<polygoncnt;i++){
			if(i != external){
				paintPolygon(i,0);
				polygonsToExtrude[pol] = i;
				pol++;
			}
		}
	}
	if(rule == 1){
		for(i=0;i<polygoncnt;i++){
			if(i != external){
				if(polygonList[i].w_num % 2 == 1 || 
				polygonList[i].w_num % 2 == -1 ){
					paintPolygon(i,0);
					polygonsToExtrude[pol] = i;
					pol++;
				}
			}
		}
	}
	if(rule == 2){
		for(i=0;i<polygoncnt;i++){
			if(i != external){
				if(polygonList[i].w_num != 0){
					paintPolygon(i,0);
					polygonsToExtrude[pol] = i;
					pol++;
				}
			}
		}
	}
	if(rule == 3){
		for(i=0;i<polygoncnt;i++){
			if(i != external){
				if((int) ((polygonList[i].w_num / (double) 2)+ 0.5) % 2 == 1 || 
					(int) ((polygonList[i].w_num / (double) 2)+ 0.5) % 2 == -1){
					paintPolygon(i,0);
					polygonsToExtrude[pol] = i;
					pol++;
				}
			}
		}
	}


}

void paintPolygon(int i,int z){ //doesn't work in Windows
    int j;
    GLdouble v[40][3] = {-1};
    tobj = gluNewTess();
    gluTessProperty(tobj, GLU_TESS_WINDING_NONZERO,currentWinding);

    gluTessCallback(tobj, GLU_TESS_VERTEX,
                   (void (CALLBACK*)()) glVertex3dv);
    gluTessCallback(tobj, GLU_TESS_BEGIN,
    		(void (CALLBACK*)()) beginCallback);
    gluTessCallback(tobj, GLU_TESS_END,
                   (void (CALLBACK*)()) endCallback);
    gluTessCallback(tobj, GLU_TESS_ERROR,
                   (void (CALLBACK*)()) errorCallback);
	if(option == 41){
    		glColor3f(line_red,line_green,line_blue);	
	}
	else{
    		glColor3f(fill_red,fill_green,fill_blue);
	}
    glShadeModel(GL_SMOOTH); 
    gluTessBeginPolygon(tobj, NULL);
        gluTessBeginContour(tobj); 
           for(j=0;j<polygonList[i].cnt;j++){
		
               v[j][0] = polygonList[i].polygonEdges[j].x1;
               v[j][1] = polygonList[i].polygonEdges[j].y1;
               v[j][2] = z;
		gluTessVertex(tobj, v[j], v[j]);
	       
           }
      	gluTessEndContour(tobj);
      gluTessEndPolygon(tobj);
}

int interEdgesNotExtrude(int x1,int y1,int x2,int y2){
	int i,j;

	for(i=0;i<pol;i++){
		for(j=0;j<polygonList[polygonsToExtrude[i]].cnt;j++){
			if(x1 == polygonList[polygonsToExtrude[i]].polygonEdges[j].x2 &&
			y1 == polygonList[polygonsToExtrude[i]].polygonEdges[j].y2 &&
			x2 == polygonList[polygonsToExtrude[i]].polygonEdges[j].x1 &&
			y2 == polygonList[polygonsToExtrude[i]].polygonEdges[j].y1 ){
				return 1;
			}
	
		}
	}
	return -1;

}

void extrude(){
	int i,j;
    glMatrixMode (GL_PROJECTION);
    glLoadIdentity(); 
    glMatrixMode (GL_MODELVIEW);
	glLoadIdentity();
    glPushMatrix();
 
    glOrtho(-600.0, 600.0, -500.0, 500.0, -3000.0, 3000.0);
    gluLookAt(alpha,250.0 , beta ,  300.0,  250.0, -100.0,  0,  1,  0);
    glColor3f(line_red,line_green,line_blue);
    glBegin(GL_LINES);
	for(i = 0;i<pol;i++){
		for(j=0;j<polygonList[polygonsToExtrude[i]].cnt;j++){
			if(interEdgesNotExtrude(polygonList[polygonsToExtrude[i]].polygonEdges[j].x1 ,polygonList[polygonsToExtrude[i]].polygonEdges[j].y1,
				polygonList[polygonsToExtrude[i]].polygonEdges[j].x2 ,polygonList[polygonsToExtrude[i]].polygonEdges[j].y2) == -1){
				glVertex3f(polygonList[polygonsToExtrude[i]].polygonEdges[j].x1 ,polygonList[polygonsToExtrude[i]].polygonEdges[j].y1 ,-50);
				glVertex3f(polygonList[polygonsToExtrude[i]].polygonEdges[j].x2 ,polygonList[polygonsToExtrude[i]].polygonEdges[j].y2 ,-50);
				glVertex3f(polygonList[polygonsToExtrude[i]].polygonEdges[j].x1 ,polygonList[polygonsToExtrude[i]].polygonEdges[j].y1 ,-100);
				glVertex3f(polygonList[polygonsToExtrude[i]].polygonEdges[j].x2 ,polygonList[polygonsToExtrude[i]].polygonEdges[j].y2 ,-100);
				glVertex3f(polygonList[polygonsToExtrude[i]].polygonEdges[j].x1 ,polygonList[polygonsToExtrude[i]].polygonEdges[j].y1 ,-50);
				glVertex3f(polygonList[polygonsToExtrude[i]].polygonEdges[j].x1 ,polygonList[polygonsToExtrude[i]].polygonEdges[j].y1 ,-100);
				glVertex3f(polygonList[polygonsToExtrude[i]].polygonEdges[j].x2 ,polygonList[polygonsToExtrude[i]].polygonEdges[j].y2 ,-50);
				glVertex3f(polygonList[polygonsToExtrude[i]].polygonEdges[j].x2 ,polygonList[polygonsToExtrude[i]].polygonEdges[j].y2 ,-100);
			}
		}
	}
	glEnd();

//	glFrontFace(GL_CCW);
	glEnable(GL_CULL_FACE);

	for(i = 0;i<pol;i++){
		for(j=0;j<polygonList[polygonsToExtrude[i]].cnt;j++){	
			if(interEdgesNotExtrude(polygonList[polygonsToExtrude[i]].polygonEdges[j].x1 ,polygonList[polygonsToExtrude[i]].polygonEdges[j].y1,
				polygonList[polygonsToExtrude[i]].polygonEdges[j].x2 ,polygonList[polygonsToExtrude[i]].polygonEdges[j].y2) == -1){
				glColor3f(fill_red,fill_green,fill_blue);
				glBegin(GL_QUADS);
				glVertex3f(polygonList[polygonsToExtrude[i]].polygonEdges[j].x1  ,polygonList[polygonsToExtrude[i]].polygonEdges[j].y1,-101);
				glVertex3f(polygonList[polygonsToExtrude[i]].polygonEdges[j].x1 ,polygonList[polygonsToExtrude[i]].polygonEdges[j].y1,-51);
				glVertex3f(polygonList[polygonsToExtrude[i]].polygonEdges[j].x2 ,polygonList[polygonsToExtrude[i]].polygonEdges[j].y2,-51);
				glVertex3f(polygonList[polygonsToExtrude[i]].polygonEdges[j].x2,polygonList[polygonsToExtrude[i]].polygonEdges[j].y2,-101);
				glEnd();
			}
		}
	}
	for(i = 0;i<pol;i++){
		for(j=0;j<polygonList[polygonsToExtrude[i]].cnt;j++){
			if(interEdgesNotExtrude(polygonList[polygonsToExtrude[i]].polygonEdges[j].x1 ,polygonList[polygonsToExtrude[i]].polygonEdges[j].y1,
				polygonList[polygonsToExtrude[i]].polygonEdges[j].x2 ,polygonList[polygonsToExtrude[i]].polygonEdges[j].y2) == -1){
				glColor3f(line_red, line_green, line_blue);
				glBegin(GL_QUADS);
				glVertex3f(polygonList[polygonsToExtrude[i]].polygonEdges[j].x1 ,polygonList[polygonsToExtrude[i]].polygonEdges[j].y1 ,-50);
				glVertex3f(polygonList[polygonsToExtrude[i]].polygonEdges[j].x1 ,polygonList[polygonsToExtrude[i]].polygonEdges[j].y1 ,-100);
				glVertex3f(polygonList[polygonsToExtrude[i]].polygonEdges[j].x2 ,polygonList[polygonsToExtrude[i]].polygonEdges[j].y2 ,-100);
				glVertex3f(polygonList[polygonsToExtrude[i]].polygonEdges[j].x2 ,polygonList[polygonsToExtrude[i]].polygonEdges[j].y2 ,-50);
				glEnd();
			}
		}
	}
	glDisable(GL_CULL_FACE);
	if(option == 41){
		for(i=0;i<pol;i++){
			paintPolygon(polygonsToExtrude[i],-50);
 			paintPolygon(polygonsToExtrude[i],-100);
			 
		}
	}	
	glPopMatrix();
}

void extrudeDifLen(){
	int i,j,increasingLength = -10;
    glMatrixMode (GL_PROJECTION);
    glLoadIdentity(); 
    glMatrixMode (GL_MODELVIEW);
	glLoadIdentity();
    glPushMatrix();
 
    glOrtho(-600.0, 600.0, -500.0, 500.0, -3000.0, 3000.0);
    gluLookAt(alpha,250.0 , beta ,  300.0,  250.0, -100.0,  0,  1,  0);
    glColor3f(line_red,line_green,line_blue);
    for(j = 0;j < pol;j++){
   	glBegin(GL_LINES);
    	for(i=0;i<polygonList[polygonsToExtrude[j]].cnt;i++){
			if(interEdgesNotExtrude(polygonList[polygonsToExtrude[j]].polygonEdges[i].x1 ,polygonList[polygonsToExtrude[j]].polygonEdges[i].y1,
				polygonList[polygonsToExtrude[j]].polygonEdges[i].x2 ,polygonList[polygonsToExtrude[j]].polygonEdges[i].y2) == -1){
				glVertex3f(polygonList[polygonsToExtrude[j]].polygonEdges[i].x1 ,polygonList[polygonsToExtrude[j]].polygonEdges[i].y1 ,increasingLength);
				glVertex3f(polygonList[polygonsToExtrude[j]].polygonEdges[i].x2 ,polygonList[polygonsToExtrude[j]].polygonEdges[i].y2 ,increasingLength);
				glVertex3f(polygonList[polygonsToExtrude[j]].polygonEdges[i].x1 ,polygonList[polygonsToExtrude[j]].polygonEdges[i].y1 ,increasingLength*4);
				glVertex3f(polygonList[polygonsToExtrude[j]].polygonEdges[i].x2 ,polygonList[polygonsToExtrude[j]].polygonEdges[i].y2 ,increasingLength*4);
				glVertex3f(polygonList[polygonsToExtrude[j]].polygonEdges[i].x1 ,polygonList[polygonsToExtrude[j]].polygonEdges[i].y1 ,increasingLength);
				glVertex3f(polygonList[polygonsToExtrude[j]].polygonEdges[i].x1 ,polygonList[polygonsToExtrude[j]].polygonEdges[i].y1 ,increasingLength*4);
				glVertex3f(polygonList[polygonsToExtrude[j]].polygonEdges[i].x2 ,polygonList[polygonsToExtrude[j]].polygonEdges[i].y2 ,increasingLength);
				glVertex3f(polygonList[polygonsToExtrude[j]].polygonEdges[i].x2 ,polygonList[polygonsToExtrude[j]].polygonEdges[i].y2 ,increasingLength*4);
			}

	}
	increasingLength -= 10;
    	glEnd();
     }
	glEnable(GL_CULL_FACE);
	increasingLength = -10;
	for(i = 0;i<pol;i++){
		for(j=0;j<polygonList[polygonsToExtrude[i]].cnt;j++){	
			if(interEdgesNotExtrude(polygonList[polygonsToExtrude[i]].polygonEdges[j].x1 ,polygonList[polygonsToExtrude[i]].polygonEdges[j].y1,
				polygonList[polygonsToExtrude[i]].polygonEdges[j].x2 ,polygonList[polygonsToExtrude[i]].polygonEdges[j].y2) == -1){
				glColor3f(fill_red,fill_green,fill_blue);
				glBegin(GL_QUADS);
				glVertex3f(polygonList[polygonsToExtrude[i]].polygonEdges[j].x1  ,polygonList[polygonsToExtrude[i]].polygonEdges[j].y1,increasingLength*4);
				glVertex3f(polygonList[polygonsToExtrude[i]].polygonEdges[j].x1 ,polygonList[polygonsToExtrude[i]].polygonEdges[j].y1,increasingLength);
				glVertex3f(polygonList[polygonsToExtrude[i]].polygonEdges[j].x2 ,polygonList[polygonsToExtrude[i]].polygonEdges[j].y2,increasingLength);
				glVertex3f(polygonList[polygonsToExtrude[i]].polygonEdges[j].x2,polygonList[polygonsToExtrude[i]].polygonEdges[j].y2,increasingLength*4);
				glEnd();
			}
		}
		increasingLength -= 10;
	}
	increasingLength = -10;
	for(i = 0;i<pol;i++){
		for(j=0;j<polygonList[polygonsToExtrude[i]].cnt;j++){
			if(interEdgesNotExtrude(polygonList[polygonsToExtrude[i]].polygonEdges[j].x1 ,polygonList[polygonsToExtrude[i]].polygonEdges[j].y1,
				polygonList[polygonsToExtrude[i]].polygonEdges[j].x2 ,polygonList[polygonsToExtrude[i]].polygonEdges[j].y2) == -1){
				glColor3f(line_red, line_green, line_blue);
				glBegin(GL_QUADS);
				glVertex3f(polygonList[polygonsToExtrude[i]].polygonEdges[j].x1 ,polygonList[polygonsToExtrude[i]].polygonEdges[j].y1 ,increasingLength);
				glVertex3f(polygonList[polygonsToExtrude[i]].polygonEdges[j].x1 ,polygonList[polygonsToExtrude[i]].polygonEdges[j].y1 ,increasingLength*4);
				glVertex3f(polygonList[polygonsToExtrude[i]].polygonEdges[j].x2 ,polygonList[polygonsToExtrude[i]].polygonEdges[j].y2 ,increasingLength*4);
				glVertex3f(polygonList[polygonsToExtrude[i]].polygonEdges[j].x2 ,polygonList[polygonsToExtrude[i]].polygonEdges[j].y2 ,increasingLength);
				glEnd();
			}
		}
		increasingLength -= 10;
	}
	glDisable(GL_CULL_FACE);

	if(option == 41){
		increasingLength = -10;
		for(i=0;i<pol;i++){
			paintPolygon(polygonsToExtrude[i],increasingLength);
 			paintPolygon(polygonsToExtrude[i],increasingLength*4);
			increasingLength -= 10;
		}
	}	
     glPopMatrix();
}

static void sfunc(int k, int x, int y)
{
	//double alpha;
// Here we can process function keys and other special key events
   switch (k)
   {

// Rotate to the left
	  case GLUT_KEY_LEFT:
		alpha = sin(angle)*600;
		beta = -cos(angle)*600;
		angle -= 0.015f;
		  display();
	  break;

// Rotate to the right
	  case GLUT_KEY_RIGHT:
		alpha = sin(angle)*600;
		beta = -cos(angle)*600;
		angle += 0.015f;
		  display();
	  break;
   }
}




void pickFillColor(){
	if(option==2){
		fill_red=1.0;
		fill_green=1.0;
		fill_blue=1.0;
	}else if(option==3){
		fill_red=0.0;
		fill_green=0.0;
		fill_blue=0.0;
	}else if(option==4){
		fill_red=1.0;
		fill_green=0.0;
		fill_blue=0.0;
	}else if(option==5){
		fill_red=0.0;
		fill_green=1.0;
		fill_blue=0.0;
	}else if(option==6){
		fill_red=0.0;
		fill_green=0.0;
		fill_blue=1.0;
	}else if(option==7){
		fill_red=1.0;
		fill_green=1.0;
		fill_blue=0.0;
	}else if(option==8){
		fill_red=1.0;
		fill_green=0.0;
		fill_blue=1.0;
	}else if(option==9){
		fill_red=0.0;
		fill_green=1.0;
		fill_blue=1.0;
	}else if(option==10){
		fill_red=1.0;
		fill_green=0.5;
		fill_blue=0.0;
	}else if(option==11){
		fill_red= 0.647059;
		fill_green=0.164706;
		fill_blue=0.164706;
	}else if(option==12){
		fill_red=0.737255;
		fill_green=0.560784;
		fill_blue=0.560784;
	}else if(option==13){
		fill_red=0.309804;
		fill_green= 0.184314;
		fill_blue=0.309804;
	}else if(option==14){
		fill_red= 0.8;
		fill_green=0.498039;
		fill_blue=0.196078;
	}else if(option==15){
		fill_red=0.90;
		fill_green=0.91;
		fill_blue=0.98;
	}else if(option==16){
		fill_red=0.30;
		fill_green=0.30;
		fill_blue=1.0;
	}else if(option==17){
		fill_red=0.6;
		fill_green= 0.8;
		fill_blue=0.196078;
	}
}


void pickLineColor(){
	if(option==18){
		line_red=1.0;
		line_green=1.0;
		line_blue=1.0;
	}else if(option==19){
		line_red=0.0;
		line_green=0.0;
		line_blue=0.0;
	}else if(option==20){
		line_red=1.0;
		line_green=0.0;
		line_blue=0.0;
	}else if(option==21){
		line_red=0.0;
		line_green=1.0;
		line_blue=0.0;
	}else if(option==22){
		line_red=0.0;
		line_green=0.0;
		line_blue=1.0;
	}else if(option==23){
		line_red=1.0;
		line_green=1.0;
		line_blue=0.0;
	}else if(option==24){
		line_red=1.0;
		line_green=0.0;
		line_blue=1.0;
	}else if(option==25){
		line_red=0.0;
		line_green=1.0;
		line_blue=1.0;
	}else if(option==26){
		line_red=1.0;
		line_green=0.5;
		line_blue=0.0;
	}else if(option==27){
		line_red= 0.647059;
		line_green=0.164706;
		line_blue=0.164706;
	}else if(option==28){
		line_red=0.737255;
		line_green=0.560784;
		line_blue=0.560784;
	}else if(option==29){
		line_red=0.309804;
		line_green= 0.184314;
		line_blue=0.309804;
	}else if(option==30){
		line_red= 0.8;
		line_green=0.498039;
		line_blue=0.196078;
	}else if(option==31){
		line_red=0.90;
		line_green=0.91;
		line_blue=0.98;
	}else if(option==32){
		line_red=0.30;
		line_green=0.30;
		line_blue=1.0;
	}else if(option==33){
		line_red=0.6;
		line_green= 0.8;
		line_blue=0.196078;
	}
}

void drawChangedPolygon(){
	int i,j;
	glColor3f(line_red, line_green, line_blue);
	for(i=0;i<(cordpnt+1);i++){
		glPointSize(POINT_SIZE);
		glBegin(GL_POINTS); //zwgrafisma tou x,y shmeioy
  		glVertex2f(cord[i][0],cord[i][1]);
		glEnd();


		if(cordpnt>0 && i>0){
			glBegin(GL_LINES);
			glVertex2f(cord[i-1][0], cord[i-1][1]);
			glVertex2f(cord[i][0], cord[i][1]);
			glEnd();

		}
	}
	//enwsh teleutaiou shmeiou me to prwto
	//sxediazw tis grammes (pernoun to xrwma pou exei epilexteis, vlepe panw)
	glBegin(GL_LINES);
	glVertex2f(cord[cordpnt][0], cord[cordpnt][1]);
	glVertex2f(cord[0][0], cord[0][1]);
	glEnd();

	for(i=0;i<cordpnt+1;i++){
		if(i==cordpnt){
			edgeList[i].x1 = cord[i][0];
	        edgeList[i].x2 = cord[0][0];
	        edgeList[i].y1 = cord[i][1];
	        edgeList[i].y2 = cord[0][1];
	        edgeList[i].cnt=0;

		}else{
			edgeList[i].x1 = cord[i][0];
	        edgeList[i].x2 = cord[i+1][0];
	        edgeList[i].y1 = cord[i][1];
	        edgeList[i].y2 = cord[i+1][1];
	        edgeList[i].cnt=0;
		}
	}

	findIntersectionPoints2();
	divideEdges();
	findPolygons();
	findWindingNumbers();
	pickPolygonsToFill();


	option=NONE;
	createMenu(); //3ana ftiaxnw to menu giati aferouse(gia kapoio logo) ta options tou xrwmmatos
//	glutAttachMenu(GLUT_RIGHT_BUTTON);
	//glPopMatrix();
	for(i=0;i<(cordpnt+1);i++){
		printf("%d--->(%d,%d)\n",i,cord[i][0], cord[i][1]);
	}
	
	printf("Polygon with w_num = -1000 is the external.\n");
	printf("**************************\n");

	for(i=0;i<polygoncnt;i++){
		printf("Polygon %d with w_num=%d\n",i,polygonList[i].w_num);
		for(j=0;j<polygonList[i].cnt;j++){
			 printf("(x%d: %d, y%d: %d) (x%d: %d, y%d: %d) --- ",j,polygonList[i].polygonEdges[j].x1,j,
                   polygonList[i].polygonEdges[j].y1,j+1,polygonList[i].polygonEdges[j].x2,j+1,polygonList[i].polygonEdges[j].y2);
		}
		printf("\n");
		printf("**************************\n");
	}
}

void drawPolygon(){
	int i,j;
	alpha = 600;	//reseting the camera angles (for extrudes)
	beta = 600;
	//glPushMatrix();
	glColor3f(line_red, line_green, line_blue);
	if(button== GLUT_LEFT_BUTTON){
		if(option==1){
			cordpnt++;
			cord[cordpnt][0]=x; //perasma tou x,y sto pinaka twn akmwn
			cord[cordpnt][1]=y;
		}

		for(i=0;i<(cordpnt+1);i++){
			glPointSize(POINT_SIZE);
			glBegin(GL_POINTS); //zwgrafisma tou x,y shmeioy
	   		glVertex2f(cord[i][0],cord[i][1]);
			glEnd();
			if(cordpnt>0 && i>0){
				glBegin(GL_LINES);
				glVertex2f(cord[i-1][0], cord[i-1][1]);
				glVertex2f(cord[i][0], cord[i][1]);
				glEnd();

			}
		}


	}else if(button== GLUT_RIGHT_BUTTON){
		for(i=0;i<(cordpnt+1);i++){
			glPointSize(10);
			glBegin(GL_POINTS); //zwgrafisma tou x,y shmeioy
	   		glVertex2f(cord[i][0],cord[i][1]);
			glEnd();


			if(cordpnt>0 && i>0){
				glBegin(GL_LINES);
				glVertex2f(cord[i-1][0], cord[i-1][1]);
				glVertex2f(cord[i][0], cord[i][1]);
				glEnd();

			}
		}
		//enwsh teleutaiou shmeiou me to prwto

		//sxediazw tis grammes (pernoun to xrwma pou exei epilexteis, vlepe panw)
		glBegin(GL_LINES);
		glVertex2f(cord[cordpnt][0], cord[cordpnt][1]);
		glVertex2f(cord[0][0], cord[0][1]);
		glEnd();

		//findIntersectionPoints();

		for(i=0;i<cordpnt+1;i++){
			if(i==cordpnt){
				edgeList[i].x1 = cord[i][0];
	            edgeList[i].x2 = cord[0][0];
	            edgeList[i].y1 = cord[i][1];
	            edgeList[i].y2 = cord[0][1];
	            edgeList[i].cnt=0;
			}else{
				edgeList[i].x1 = cord[i][0];
	            edgeList[i].x2 = cord[i+1][0];
	            edgeList[i].y1 = cord[i][1];
	            edgeList[i].y2 = cord[i+1][1];
	            edgeList[i].cnt=0;
			}
		}
		findIntersectionPoints2();
		divideEdges();
		findPolygons();
		findWindingNumbers();
		pickPolygonsToFill();
		for(i=0;i<(cordpnt+1);i++){
			printf("%d--->(%d,%d)\n",i,cord[i][0], cord[i][1]);
		}
	
		printf("Polygon with w_num = -1000 is the external.\n");
		printf("**************************\n");

		for(i=0;i<polygoncnt;i++){
			printf("Polygon %d with w_num=%d\n",i,polygonList[i].w_num);
			for(j=0;j<polygonList[i].cnt;j++){
				 printf("(x%d: %d, y%d: %d) (x%d: %d, y%d: %d) --- ",j,polygonList[i].polygonEdges[j].x1,j,
		           polygonList[i].polygonEdges[j].y1,j+1,polygonList[i].polygonEdges[j].x2,j+1,polygonList[i].polygonEdges[j].y2);
			}
			printf("\n");
			printf("**************************\n");
		}
		
			//gia na upostirizetai to de3i click-->menu oso den eimai sthn Polygon
			option=NONE;
			createMenu(); 
	}
	


}

int MovePoint(){
	int i;
	for(i=0;i<(cordpnt+1);i++){
		if(xi>=(cord[i][0]-POINT_SIZE) && xi<=(cord[i][0]+POINT_SIZE) && yi>=(cord[i][1]-POINT_SIZE) && yi<=(cord[i][1]+POINT_SIZE)){
			return i;
		}
	}
	return -1;

}

void menu(int value){
	if(value == 0){  //exit
		glutDestroyWindow(window);
		exit(0);
	}else if(value == 1){ //polygon
		Init();
		cordpnt=-1;
		option=value;
	}else if(value>=2 && value<=17){ //fill color
		option=value;
		pickFillColor();
	}else if(value>=18 && value<=33){ //line color
		option=value;
		pickLineColor();
	}else if(value==34){ //move
		option=value;
	}else if(value==35){ //extrude
		extrudeType = 0;
		option=value;
		display();
	}
	else if(value==36){ //extrude w different length
		extrudeType = 1;
		option=value;
		display();
	}
	else if(value==37){   //RULES
		option = value;
		rule = 0;
	}
	else if(value==38){ 
		option=value;
		rule = 1;
	}
	else if(value==39){ 
		option=value;
		rule = 2;
	}
	else if(value==40){ 
		option=value;
		rule = 3;
	}else if(value==41){ //Close Polygon
		option=value;
		display();	
	}

}

void createMenu(void){
	int action,polygon,fillColor,lineColor,rules;

	action = glutCreateMenu(menu);
	glutAddMenuEntry("EXIT",0);
	glutAddMenuEntry("Polygon",1);
	glutAddMenuEntry("Move",34);
	glutAddMenuEntry("EXTRUDE",35);
	glutAddMenuEntry("EXTR PLUS",36);
	glutAddMenuEntry("Close Polygon",41);
	glutCreateMenu(menu);

	fillColor=glutCreateMenu(menu);
	glutAddMenuEntry("WHITE",2);
	glutAddMenuEntry("BLACK",3);
	glutAddMenuEntry("RED",4);
	glutAddMenuEntry("GREEN",5);
	glutAddMenuEntry("BLUE",6);
	glutAddMenuEntry("YELLOW",7);
	glutAddMenuEntry("PURPLE",8);
	glutAddMenuEntry("CYAN",9);
	glutAddMenuEntry("ORANGE",10);
	glutAddMenuEntry("BROWN",11);
	glutAddMenuEntry("PINK",12);
	glutAddMenuEntry("VIOLET",13);
	glutAddMenuEntry("GOLD",14);
	glutAddMenuEntry("SILVER",15);
	glutAddMenuEntry("NEON_BLUE",16);
	glutAddMenuEntry("YELLOW_GREEN",17);
	//glutCreateMenu(menu);

	lineColor=glutCreateMenu(menu);
	glutAddMenuEntry("WHITE",18);
	glutAddMenuEntry("BLACK",19);
	glutAddMenuEntry("RED",20);
	glutAddMenuEntry("GREEN",21);
	glutAddMenuEntry("BLUE",22);
	glutAddMenuEntry("YELLOW",23);
	glutAddMenuEntry("PURPLE",24);
	glutAddMenuEntry("CYAN",25);
	glutAddMenuEntry("ORANGE",26);
	glutAddMenuEntry("BROWN",27);
	glutAddMenuEntry("PINK",28);
	glutAddMenuEntry("VIOLET",29);
	glutAddMenuEntry("GOLD",30);
	glutAddMenuEntry("SILVER",31);
	glutAddMenuEntry("NEON_BLUE",32);
	glutAddMenuEntry("YELLOW_GREEN",33);
	glutCreateMenu(menu);

	rules = glutCreateMenu(menu);
	glutAddMenuEntry("Rule 0",37);
	glutAddMenuEntry("Rule 1",38);
	glutAddMenuEntry("Rule 2",39);
	glutAddMenuEntry("Rule 3",40);
	glutCreateMenu(menu);

	glutAddSubMenu("ACTION",action);
	glutAddSubMenu("FILL_COLOR",fillColor);
	glutAddSubMenu("LINE_COLOR",lineColor);
	glutAddSubMenu("Rules",rules);
	glutAttachMenu(GLUT_RIGHT_BUTTON);

}

void windowMouse(int btn, int state, int a, int b)
{
	if (btn == GLUT_LEFT_BUTTON){

		button=btn;
  	}else if(btn == GLUT_RIGHT_BUTTON){
  		button=btn;
  	}


 	if (state == GLUT_DOWN){ //Mouse DOWN
  		if(option==34){
  			xi=a;
  			yi=500-b;

  			found=MovePoint();
			if(found==-1){
				createMenu();
				option = NONE;
			}
  		}
		if(option==1){
			glutDetachMenu(GLUT_RIGHT_BUTTON);

			x=a;
			y=500-b;
			display();
		}
	}else{ //mouse UP
		if(option==34){
			glutDetachMenu(GLUT_RIGHT_BUTTON);
			if(found>-1){
				cord[found][0]=a;
				cord[found][1]=500-b;
				display();
				found = -1;

				
			}
		}
	}
}

void display(){

	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);	// Clear Screen

	if(option==1){
		drawPolygon();
		glutSwapBuffers();
	}else if(option==34 && found!=-1){
		drawChangedPolygon();
		glutSwapBuffers();
	}else if(option==35){
		extrude();
		glutSwapBuffers();
	}else if(option==36){
		extrudeDifLen();
		glutSwapBuffers();
	}else if(option==41){
		if(extrudeType == 0){
			extrude();
			glutSwapBuffers();
		}
		else if(extrudeType == 1){
			extrudeDifLen();
			glutSwapBuffers();
		}
	}else if(firstTimeOpened == 0){
		glutSwapBuffers();
		firstTimeOpened = 1;
	}
}


void Init(void)
{	
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
  	glClearColor(1, 1, 1, 1);	//dinei sto parathuro aspro,telautaia metavlhth dhlwnei diafaneia
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	gluOrtho2D(0.0, 600.0,0.0, 500.0);
	glLineWidth(2);
	glutSwapBuffers();
}

int main( int argc, char **argv)
{
   	glutInit(&argc, argv);
	glutInitDisplayMode(GLUT_RGB | GLUT_DOUBLE | GLUT_DEPTH );
  	glutInitWindowSize (600, 500);
   	glutInitWindowPosition (100, 100);
   	window = glutCreateWindow ("MyPolygon");
   	createMenu();


   	Init();
   	glutDisplayFunc(display);
   	glutMouseFunc(windowMouse);
    glutSpecialFunc(sfunc);
   	glutMainLoop();
   	return 0;
}

