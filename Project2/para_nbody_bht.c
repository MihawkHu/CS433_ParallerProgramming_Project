// Compile with:
//
//
// To specify the number of bodies in the world, the program optionally accepts
// an integer as its first command line argument.

#include <time.h>
#include <sys/times.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <X11/Xlib.h>
#include <unistd.h>

#ifdef _OPENMP
#include <omp.h>
#endif

#define WIDTH 1024
#define HEIGHT 768
#define THREAD_NUM 4
// default number of bodies
#define DEF_NUM_BODIES 200
// gravitational constant
#define GRAV 10.0
// initial velocities are scaled by this value
#define V_SCALAR 20.0
// initial masses are scaled by this value
#define M_SCALAR 5.0
// radius scalar
#define R_SCALAR 3
// coefficient of restitution determines the elasticity of a collision: C_REST = [0,1]
//  if C_REST = 0 -> perfectly inelastic (particles stick together)
//  if C_REST = 1 -> perfectly elastic (no loss of speed)
#define C_REST 0.5
// set the iteration times
#define iteration_times 200
// Must set 0 if run on Pi
#define NOT_RUN_ON_PI 1

struct body {
    double x, y; // position
    double vx, vy; // velocity
    double m; // mass
    double r; // radius of the particle
};

struct world {
    struct body *bodies;
    int num_bodies;
};

struct node{
	int level;
	double mass_x;
	double mass_y;
    double xmin,xmax,ymin,ymax;
    double totalMass;
	int sons_num;
	int isLeaf;
    int *sons;
	struct node *Son[4];
};
int flag = 1;
double ALPHA = 0.000001;
clock_t total_time = 0;
//total_time.sec = 0;
//total_time.usec = 0;

struct qNode {
    struct qNode *next;
    struct node *nn;
};

struct queue{
    struct qNode *start;
    struct qNode *end;
    int size;
};

void queueInit(struct queue **q) {
    (*q) = (struct queue*)malloc(sizeof(struct queue));
    (*q)->start = (*q)->end = NULL;
    (*q)->size = 0;
}

void push(struct queue **q, struct node **n) {
  	if ((*q)->size == 0) {
        (*q)->start = (struct qNode*)malloc(sizeof(struct qNode));
        (*q)->start->nn = (*n);
        (*q)->start->next = NULL;
        (*q)->end = (*q)->start;
    }
    else {
        struct qNode *temp = (struct qNode*)malloc(sizeof(struct qNode));
        temp->next = NULL;
        temp->nn = (*n);
        (*q)->end->next = temp;
        (*q)->end = temp;
    }
    (*q)->size++;
}

struct node* pop(struct queue **q) {
    if ((*q)->size == 0) {
        printf("Nothing in the queue!\n");
        return NULL;
    }

    struct qNode *temp = (*q)->start;
    struct node *res = temp->nn;
    (*q)->start = (*q)->start->next;

    free(temp);
    (*q)->size--;
    return res;
}

int isEmpty(struct queue **q) {
    return (*q)->size;
}

/* This function initializes a new BHT*/
void BHT_construct_root(double x, double y, struct world *world, struct node *Root){
    int i;

    Root->xmax = x;
    Root->ymax = y;
    Root->xmin = Root->ymin = 0;
	Root->level = 0;
	Root->sons_num = world->num_bodies;
    Root->Son[0] = NULL;
    Root->Son[1] = NULL;
    Root->Son[2] = NULL;
    Root->Son[3] = NULL;
    Root->isLeaf = 0;
    Root->sons = malloc(sizeof(int) * Root->sons_num);
    for(i = 0; i < world->num_bodies; i++){
        Root->sons[i] = i;
    }
}

void node_construct_2(struct node *tempNode, struct node *son){
    son->level = tempNode->level + 1;
    son->isLeaf = 1;
    son->sons_num = 0;
    son->Son[0] = son->Son[1] = son->Son[2] = son->Son[3] = NULL;
}

void node_construct_1(struct node *tempNode){
    int i;
    for(i = 0; i < 4; i++){
        tempNode->Son[i]->xmin = tempNode->xmin + (tempNode->xmax - tempNode->xmin) * (i%2) / 2;
        tempNode->Son[i]->xmax = tempNode->xmax - (tempNode->xmax - tempNode->xmin) * ((i+1)%2) / 2;
        tempNode->Son[i]->ymin = tempNode->ymin + (tempNode->ymax - tempNode->ymin) * (i%2) / 2;
        tempNode->Son[i]->ymax = tempNode->ymax - (tempNode->ymax - tempNode->ymin) * ((i+1)%2) / 2;
        node_construct_2(tempNode, tempNode->Son[i]);
    }
}

void BHT_construct(struct world *world,struct node *tempNode){
    int i, j;
    for(i = 0 ; i < 4 ;i ++) {
        tempNode->Son[i] = malloc(sizeof(struct node));
        tempNode->Son[i]->sons = malloc(sizeof(int)*tempNode->sons_num);
    }
    node_construct_1(tempNode);
    for(i = 0; i < tempNode->sons_num; i++){
        j = tempNode->sons[i];
        if((world->bodies[j].x > tempNode->xmin) & (world->bodies[j].x <= (tempNode->xmin + tempNode->xmax)/2)){
            if((world->bodies[j].y > (tempNode->ymin + tempNode->ymax)/2) & (world->bodies[j].y <= tempNode->ymax))
                tempNode->Son[0]->sons[tempNode->Son[0]->sons_num++] = j;
            else
                tempNode->Son[2]->sons[tempNode->Son[2]->sons_num++] = j;
        }
        else{
            if((world->bodies[j].y > (tempNode->ymin + tempNode->ymax)/2) & (world->bodies[j].y <= tempNode->ymax))
                tempNode->Son[1]->sons[tempNode->Son[1]->sons_num++] = j;
            else
                tempNode->Son[3]->sons[tempNode->Son[3]->sons_num++] = j;
        }
    }
    for(i = 0; i < 4; i++){
        if(tempNode->Son[i]->sons_num == 0){
            free(tempNode->Son[i]);
            tempNode->Son[i] = NULL;
        }
        else if(tempNode->Son[i]->sons_num > 1){
            tempNode->Son[i]->isLeaf = 0;
            BHT_construct(world, tempNode->Son[i]);
        }
    }
}

void BHT_complete(struct world *world,struct node *tempNode) {
    int i;
    double tempx, tempy;
    if(tempNode->isLeaf == 1) {
        tempNode->totalMass = world->bodies[tempNode->sons[0]].m;
        tempNode->mass_x = world->bodies[tempNode->sons[0]].x;
        tempNode->mass_y = world->bodies[tempNode->sons[0]].y;
        return;
    }
    tempNode->totalMass = 0;
    for(i = 0; i < 4 ; i++) {
        if(tempNode->Son[i] != NULL) {
            BHT_complete(world, tempNode->Son[i]);
            tempNode->totalMass += tempNode->Son[i]->totalMass;
        }
    }
    for(i = 0; i < 4 ; i++) {
        if(tempNode->Son[i] != NULL) {
            tempx += (tempNode->Son[i]->mass_x * tempNode->Son[i]->totalMass) / tempNode->totalMass;
            tempy += (tempNode->Son[i]->mass_y * tempNode->Son[i]->totalMass) / tempNode->totalMass;
        }
    }
    tempNode->mass_x = tempx;
    tempNode->mass_y = tempy;
}

int isContainedBy(int x, int *arr, int length) {
    int i;
    for (i = 0; i < length; i++) {
        if(arr[i] == x)
            return 1;
    }
    return 0;
}

void BHT_Force_node(int num, struct world *world, struct node *Root, double *force_x, double *force_y, struct queue *myQueue) {
    int i;
    double distance;
    double diff_x, diff_y, d_cubed;
    struct node *tmp;
    while(isEmpty(&myQueue) != 0) {
        tmp = pop(&myQueue);

        if(isContainedBy(num, tmp->sons, tmp->sons_num) == 1) {
            for(i = 0; i < 4 ; i++) {
                if(tmp->Son[i] != NULL)
                    push(&myQueue, &(tmp->Son[i]));
            }
            continue;
        }

        if(tmp->isLeaf == 1) {
            diff_x = tmp->mass_x - world->bodies[num].x;
            diff_y = tmp->mass_y - world->bodies[num].y;
            distance = sqrt(diff_x*diff_x + diff_y*diff_y);
            if(distance < 25) distance = 25;
            d_cubed = distance*distance*distance;

            force_x[num] += GRAV * (world->bodies[num].m * tmp->totalMass
                   / d_cubed) * diff_x;
            force_y[num] += GRAV * (world->bodies[num].m * tmp->totalMass
                   / d_cubed) * diff_y;
            continue;
        }

        diff_x = tmp->mass_x - world->bodies[num].x;
        diff_y = tmp->mass_y - world->bodies[num].y;
        distance = sqrt(diff_x*diff_x + diff_y*diff_y);
        if(distance < 25) distance = 25;
        d_cubed = distance*distance*distance;
        if((WIDTH / (tmp->level + 1)) / distance  < ALPHA) {
            force_x[num] += GRAV * (world->bodies[num].m * tmp->totalMass
                   / d_cubed) * diff_x;
            force_y[num] += GRAV * (world->bodies[num].m * tmp->totalMass
                   / d_cubed) * diff_y;
        }
        else {
            for(i = 0; i < 4 ; i++) {
                if(tmp->Son[i] != NULL)
                    push(&myQueue, &(tmp->Son[i]));
            }
        }
    }
    return;
}

void BHT_Force(struct world *world, struct node *Root, double *force_x, double *force_y) {
    int i;
    struct queue *myQueue;
    queueInit(&myQueue);
#   pragma omp for
    for(i = 0; i < world->num_bodies; i++) {
        push(&myQueue, &Root);
        BHT_Force_node(i, world, Root, force_x, force_y, myQueue);
    }
}

void BHT_deconstruct(struct node *Root) {
    int i;
    if(Root != NULL)
        for(i = 0; i < 4 ; i++) {
            if(Root->Son[i] != NULL)
                BHT_deconstruct(Root->Son[i]);
        }
    if(Root->sons != NULL)
        free(Root->sons);
    if(Root != NULL)
        free(Root);
}

void position_step_para(struct world *world, double time_res) {
    int i, j;
    double d, d_cubed, diff_x, diff_y;
    /* The forces array stores the x and y components of the total force acting
     * on each body. The forces are index like this:
     *     F on body i in the x dir = F_x[i]
     *     F on body i in the y dir = F_y[i] */
    double *force_x = (double*)malloc(sizeof(double) * world->num_bodies);
    double *force_y = (double*)malloc(sizeof(double) * world->num_bodies);
    // initialize all forces to zero
    force_x = memset(force_x, 0, sizeof(double) * world->num_bodies);
    force_y = memset(force_y, 0, sizeof(double) * world->num_bodies);

     struct node *myBHT;
     myBHT = malloc(sizeof(struct node));
     BHT_construct_root(1024, 768, world, myBHT);
     BHT_construct(world, myBHT);
     BHT_complete(world, myBHT);

     BHT_Force(world, myBHT, force_x, force_y);
     BHT_deconstruct(myBHT);

     for (i = 0; i < world->num_bodies; i++) {
         // Update velocities
         world->bodies[i].vx += force_x[i] * time_res / world->bodies[i].m;
         world->bodies[i].vy += force_y[i] * time_res / world->bodies[i].m;

         // Update positions
         world->bodies[i].x += world->bodies[i].vx * time_res;
         world->bodies[i].y += world->bodies[i].vy * time_res;
     }
}

/* This function initializes each particle's mass, velocity and position */
struct world* create_world(int num_bodies) {
    struct world *world = malloc(sizeof(struct world));

    world->num_bodies = num_bodies;
    world->bodies = malloc(sizeof(struct body)*num_bodies);

    int i = 0;
    double x;
    double y;
    double rc;

    int min_dim = (WIDTH < HEIGHT) ? WIDTH : HEIGHT;

    while (i<num_bodies) {
        x = drand48() * WIDTH;
        y = drand48() * HEIGHT;
        rc = sqrt((WIDTH/2-x)*(WIDTH/2-x) + (y-HEIGHT/2)*(y-HEIGHT/2));
        if (rc <= min_dim/2) {
            world->bodies[i].x = x;
            world->bodies[i].y = y;

            world->bodies[i].vx = V_SCALAR * (y-HEIGHT/2) / rc;
            world->bodies[i].vy = V_SCALAR * (WIDTH/2-x) / rc;

            world->bodies[i].m = (1 / (0.025 + drand48())) * M_SCALAR;
            world->bodies[i].r = sqrt(world->bodies[i].m / M_PI) * R_SCALAR;
            i++;
        }
    }
    return world;
}

// set the foreground color given RGB values between 0..255.
void set_color(Display *disp, GC gc, int r, int g, int b){
  unsigned long int p ;

  if (r < 0) r = 0; else if (r > 255) r = 255;
  if (g < 0) g = 0; else if (g > 255) g = 255;
  if (b < 0) b = 0; else if (b > 255) b = 255;

  p = (r << 16) | (g  << 8) | (b) ;

  XSetForeground(disp, gc, p) ;
}


/* This function updates the screen with the new positions of each particle */
void draw_world(Display *disp, Pixmap back_buf, GC gc, struct world *world) {
    int i;
    double x, y, r, r2;

    // we turn off aliasing for faster draws
    set_color(disp, gc, 255, 255, 255);
    XFillRectangle(disp, back_buf, gc, 0, 0, WIDTH, HEIGHT);

    for (i = 0; i < world->num_bodies; i++) {
        r = world->bodies[i].r;
        x = world->bodies[i].x - r;
        y = world->bodies[i].y - r;
        r2 = r + r;

        // draw body
        set_color(disp, gc, 255*7/10, 255*7/10, 255*7/10);
        XFillArc(disp, back_buf, gc, x, y, r2, r2, 0, 360*64);
        set_color(disp, gc, 0, 0, 0);
        XDrawArc(disp, back_buf, gc, x, y, r2, r2, 0, 360*64);
    }
}

void collision_step(struct world *world) {
    int a, b;
    double r, x, y, vx, vy;

    // Impose screen boundaries by reversing direction if body is off screen
    for (a = 0; a < world->num_bodies; a++) {
        r = world->bodies[a].r;
        x = world->bodies[a].x;
        y = world->bodies[a].y;
        vx = world->bodies[a].vx;
        vy = world->bodies[a].vy;

        if (x-r < 0) { // left edge
            if (vx < 0) { world->bodies[a].vx = -C_REST * vx; }
            world->bodies[a].x = r;
        } else if (x+r > WIDTH) { // right edge
            if (vx > 0) { world->bodies[a].vx = -C_REST * vx; }
            world->bodies[a].x = WIDTH - r;
        }

        if (y-r < 0) { // bottom edge
            if (vy < 0) { world->bodies[a].vy = -C_REST * vy; }
            world->bodies[a].y = r;
        } else if (y+r > HEIGHT) { // top edge
            if (vy > 0) { world->bodies[a].vy = -C_REST * vy; }
            world->bodies[a].y = HEIGHT - r;
        }
    }
}

void step_world(struct world *world, double time_res) {

	struct tms ttt;
	clock_t start, end;
	start = times(&ttt);
#	pragma omp parallel num_threads(4)
    position_step_para(world, time_res);
	end = times(&ttt);
	total_time += end - start;

    collision_step(world);
}


/* Main method runs initialize() and update() */
int main(int argc, char **argv) {
	//total_time.tv_sec = 0;
	//total_time.tv_usec = 0;
    /* get num bodies from the command line */
    int num_bodies;
    num_bodies = (argc == 2) ? atoi(argv[1]) : DEF_NUM_BODIES;
    printf("Universe has %d bodies.\n", num_bodies);

    /* set up the universe */
    time_t cur_time;
    time(&cur_time);
    srand48((long)cur_time); // seed the RNG used in create_world
    struct world *world = create_world(num_bodies);

    /* set up graphics using Xlib */
#if NOT_RUN_ON_PI
    Display *disp = XOpenDisplay(NULL);
    int scr = DefaultScreen(disp);
    Window win = XCreateSimpleWindow(
            disp,
            RootWindow(disp, scr),
            0, 0,
            WIDTH, HEIGHT,
            0,
            BlackPixel(disp, scr), WhitePixel(disp, scr));
    XStoreName(disp, win, "N-Body Simulator");

    Pixmap back_buf = XCreatePixmap(disp, RootWindow(disp, scr),
            WIDTH, HEIGHT, DefaultDepth(disp, scr));
    GC gc = XCreateGC(disp, back_buf, 0, 0);

    // Make sure we're only looking for messages about closing the window
    Atom del_window = XInternAtom(disp, "WM_DELETE_WINDOW", 0);
    XSetWMProtocols(disp, win, &del_window, 1);

    XSelectInput(disp, win, StructureNotifyMask);
    XMapWindow(disp, win);
    XEvent event;
    // wait until window is mapped
    while (1) {
        XNextEvent(disp, &event);
        if (event.type == MapNotify) {
            break;
        }
    }
#endif

    struct timespec delay={0, 1000000000 / 60}; // for 60 FPS
    struct timespec remaining;
	double delta_t = 0.1;
	int ii;

    for(ii = 0; ii < iteration_times; ii++){
        // check if the window has been closed
#if NOT_RUN_ON_PI
        if (XCheckTypedEvent(disp, ClientMessage, &event)) {
            break;
        }

        // we first draw to the back buffer then copy it to the front (`win`)
        draw_world(disp, back_buf, gc, world);
        XCopyArea(disp, back_buf, win, gc, 0, 0, WIDTH, HEIGHT, 0, 0);
#endif

        step_world(world, delta_t);

		//if you want to watch the process in 60 FPS
		//nanosleep(&delay, &remaining);
    }

//	printf("Total Time = %f\n", (double)total_time.tv_sec + (double)total_time.tv_usec/1000000);
	printf("Nbody Position Calculation Time = :%lf s\n",(double)total_time / (sysconf(_SC_CLK_TCK)));

#if NOT_RUN_ON_PI
    XFreeGC(disp, gc);
    XFreePixmap(disp, back_buf);
    XDestroyWindow(disp, win);
    XCloseDisplay(disp);
#endif

    return 0;
}
