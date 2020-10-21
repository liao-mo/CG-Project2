/************************************************************************
     File:        Maze.cpp

     Author:     
                  Stephen Chenney, schenney@cs.wisc.edu
     Modifier
                  Yu-Chi Lai, yu-chi@cs.wisc.edu

     Comment:    
						(c) 2001-2002 Stephen Chenney, University of Wisconsin at Madison

						Class header file for Maze class. Manages the maze.
		

     Platform:    Visio Studio.Net 2003 (converted to 2005)

*************************************************************************/

#include "Maze.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include <FL/Fl.h>
#include <FL/fl_draw.h>
#include <Fl/gl.h>
#include <GL/glu.h>
#include <iostream>
#include "LineSeg.h"

using namespace std;

const char Maze::X = 0;
const char Maze::Y = 1;
const char Maze::Z = 2;

const float Maze::BUFFER = 0.1f;


//**********************************************************************
//
// * Constructor for the maze exception
//======================================================================
MazeException::
MazeException(const char *m)
//======================================================================
{
	message = new char[strlen(m) + 4];
	strcpy(message, m);
}


//**********************************************************************
//
// * Constructor to create the default maze
//======================================================================
Maze::
Maze(const int nx, const int ny, const float sx, const float sy)
//======================================================================
{
	// Build the connectivity structure.
	Build_Connectivity(nx, ny, sx, sy);

	// Make edges transparent to create a maze.
	Build_Maze();

	// Set the extents of the maze
	Set_Extents();

	// Default values for the viewer.
	viewer_posn[X] = viewer_posn[Y] = viewer_posn[Z] = 0.0;
	viewer_dir = 0.0;
	viewer_fov = 45.0;

	// Always start on the 0th frame.
	frame_num = 0;
}


//**********************************************************************
//
// * Construtor to read in precreated maze
//======================================================================
Maze::
Maze(const char *filename)
//======================================================================
{
	char    err_string[128];
	FILE    *f;
	int	    i;

	// Open the file
	if ( ! ( f = fopen(filename, "r") ) )
		throw new MazeException("Maze: Couldn't open file");

	// Get the total number of vertices
	if ( fscanf(f, "%d", &num_vertices) != 1 )
		throw new MazeException("Maze: Couldn't read number of vertices");

	// Read in each vertices
	vertices = new Vertex*[num_vertices];
	for ( i = 0 ; i < num_vertices ; i++ ) {
		float x, y;
		if ( fscanf(f, "%g %g", &x, &y) != 2 )	{
			sprintf(err_string, "Maze: Couldn't read vertex number %d", i);
			throw new MazeException(err_string);
		}
		vertices[i] = new Vertex(i, x, y);
	}

	// Get the number of edges
	if ( fscanf(f, "%d", &num_edges) != 1 )
		throw new MazeException("Maze: Couldn't read number of edges");

	// read in all edges
	edges = new Edge*[num_edges];
	for ( i = 0 ; i < num_edges ; i++ ){
		int     vs, ve, cl, cr, o;
		float	r, g, b;
		if ( fscanf(f, "%d %d %d %d %d %g %g %g",
						&vs, &ve, &cl, &cr, &o, &r, &g, &b) != 8) {
			sprintf(err_string, "Maze: Couldn't read edge number %d", i);
			throw new MazeException(err_string);
		}
		edges[i] = new Edge(i, vertices[vs], vertices[ve], r, g, b);
		edges[i]->Add_Cell((Cell*)cl, Edge::LEFT);
		edges[i]->Add_Cell((Cell*)cr, Edge::RIGHT);
		edges[i]->opaque = o ? true : false;
	}

	// Read in the number of cells
	if ( fscanf(f, "%d", &num_cells) != 1 )
		throw new MazeException("Maze: Couldn't read number of cells");


	// Read in all cells
	cells = new Cell*[num_cells];
	for ( i = 0 ; i < num_cells ; i++ )	{
		int epx, epy, emx, emy;
		if ( fscanf(f, "%d %d %d %d", &epx, &epy, &emx, &emy) != 4 ){
			sprintf(err_string, "Maze: Couldn't read cell number %d", i);
			throw new MazeException(err_string);
		}
		cells[i] = new Cell(i, epx >= 0 ? edges[epx] : NULL,
									epy >= 0 ? edges[epy] : NULL,
									emx >= 0 ? edges[emx] : NULL,
									emy >= 0 ? edges[emy] : NULL);
		if ( cells[i]->edges[0] ) {
			if ( cells[i]->edges[0]->neighbors[0] == (Cell*)i )
				cells[i]->edges[0]->neighbors[0] = cells[i];
			else if ( cells[i]->edges[0]->neighbors[1] == (Cell*)i )
				cells[i]->edges[0]->neighbors[1] = cells[i];
			else	{
				sprintf(err_string,
						  "Maze: Cell %d not one of edge %d's neighbors",
							i, cells[i]->edges[0]->index);
				throw new MazeException(err_string);
			}
		}

		if ( cells[i]->edges[1] )	{
			if ( cells[i]->edges[1]->neighbors[0] == (Cell*)i )
				cells[i]->edges[1]->neighbors[0] = cells[i];
			else if ( cells[i]->edges[1]->neighbors[1] == (Cell*)i )
				cells[i]->edges[1]->neighbors[1] = cells[i];
			else {
				sprintf(err_string,
							"Maze: Cell %d not one of edge %d's neighbors",
							i, cells[i]->edges[1]->index);
				throw new MazeException(err_string);
			}
		}
		if ( cells[i]->edges[2] ) {
			if ( cells[i]->edges[2]->neighbors[0] == (Cell*)i )
				cells[i]->edges[2]->neighbors[0] = cells[i];
			else if ( cells[i]->edges[2]->neighbors[1] == (Cell*)i )
				cells[i]->edges[2]->neighbors[1] = cells[i];
			else	{
				sprintf(err_string,
							"Maze: Cell %d not one of edge %d's neighbors",
							i, cells[i]->edges[2]->index);
				throw new MazeException(err_string);
			}
		}
		if ( cells[i]->edges[3] ) {
			if ( cells[i]->edges[3]->neighbors[0] == (Cell*)i )
				cells[i]->edges[3]->neighbors[0] = cells[i];
			else if ( cells[i]->edges[3]->neighbors[1] == (Cell*)i )
				cells[i]->edges[3]->neighbors[1] = cells[i];
			else	{
				sprintf(err_string,
							"Maze: Cell %d not one of edge %d's neighbors",
							i, cells[i]->edges[3]->index);
				throw new MazeException(err_string);
			}
		}
	}

	if ( fscanf(f, "%g %g %g %g %g",
					 &(viewer_posn[X]), &(viewer_posn[Y]), &(viewer_posn[Z]),
					 &(viewer_dir), &(viewer_fov)) != 5 )
		throw new MazeException("Maze: Error reading view information.");

	// Some edges have no neighbor on one side, so be sure to set their
	// pointers to NULL. (They were set at -1 by the save/load process.)
	for ( i = 0 ; i < num_edges ; i++ )	{
		if ( edges[i]->neighbors[0] == (Cell*)-1 )
			edges[i]->neighbors[0] = NULL;
		if ( edges[i]->neighbors[1] == (Cell*)-1 )
			edges[i]->neighbors[1] = NULL;
	}

	fclose(f);

	Set_Extents();

	// Figure out which cell the viewer is in, starting off by guessing the
	// 0th cell.
	Find_View_Cell(cells[0]);

	frame_num = 0;
}


//**********************************************************************
//
// * Destructor must free all the memory allocated.
//======================================================================
Maze::
~Maze(void)
//======================================================================
{
	int i;

	for ( i = 0 ; i < num_vertices ; i++ )
		delete vertices[i];
	delete[] vertices;

	for ( i = 0 ; i < num_edges ; i++ )
		delete edges[i];
	delete[] edges;

	for ( i = 0 ; i < num_cells ; i++ )
		delete cells[i];
	delete[] cells;
}


//**********************************************************************
//
// * Randomly generate the edge's opaque and transparency for an empty maze
//======================================================================
void Maze::
Build_Connectivity(const int num_x, const int num_y,
                   const float sx, const float sy)
//======================================================================
{
	int	i, j, k;
	int edge_i;

	// Ugly code to allocate all the memory for a new maze and to associate
	// edges with vertices and faces with edges.

	// Allocate and position the vertices.
	num_vertices = ( num_x + 1 ) * ( num_y + 1 );
	vertices = new Vertex*[num_vertices];
	k = 0;
	for ( i = 0 ; i < num_y + 1 ; i++ ) {
		for ( j = 0 ; j < num_x + 1 ; j++ )	{
			vertices[k] = new Vertex(k, j * sx, i * sy);
			k++;
		}
	}

	// Allocate the edges, and associate them with their vertices.
	// Edges in the x direction get the first num_x * ( num_y + 1 ) indices,
	// edges in the y direction get the rest.
	num_edges = (num_x+1)*num_y + (num_y+1)*num_x;
	edges = new Edge*[num_edges];
	k = 0;
	for ( i = 0 ; i < num_y + 1 ; i++ ) {
		int row = i * ( num_x + 1 );
		for ( j = 0 ; j < num_x ; j++ ) {
			int vs = row + j;
			int ve = row + j + 1;
			edges[k] = new Edge(
				k, 
				vertices[vs], 
				vertices[ve],
				rand() / (float)RAND_MAX * 0.5f + 0.25f,
				rand() / (float)RAND_MAX * 0.5f + 0.25f,
				rand() / (float)RAND_MAX * 0.5f + 0.25f
			);
			k++;
		}
	}

	edge_i = k;
	for ( i = 0 ; i < num_y ; i++ ) {
		int row = i * ( num_x + 1 );
		for ( j = 0 ; j < num_x + 1 ; j++ )	{
			int vs = row + j;
			int ve = row + j + num_x + 1;
			edges[k] = new Edge(k, vertices[vs], vertices[ve],
			rand() / (float)RAND_MAX * 0.5f + 0.25f,
			rand() / (float)RAND_MAX * 0.5f + 0.25f,
			rand() / (float)RAND_MAX * 0.5f + 0.25f);
			k++;
		}
	}

	// Allocate the cells and associate them with their edges.
	num_cells = num_x * num_y;
	cells = new Cell*[num_cells];
	k = 0;
	for ( i = 0 ; i < num_y ; i++ ) {
		int row_x = i * ( num_x + 1 );
		int row_y = i * num_x;
		for ( j = 0 ; j < num_x ; j++ )	{
			int px = edge_i + row_x + 1 + j;
			int py = row_y + j + num_x;
			int mx = edge_i + row_x + j;
			int my = row_y + j;
			cells[k] = new Cell(k, edges[px], edges[py], edges[mx], edges[my]);
			edges[px]->Add_Cell(cells[k], Edge::LEFT);
			edges[py]->Add_Cell(cells[k], Edge::RIGHT);
			edges[mx]->Add_Cell(cells[k], Edge::RIGHT);
			edges[my]->Add_Cell(cells[k], Edge::LEFT);
			k++;
		}
	}
}


//**********************************************************************
//
// * Add edges from cell to the set that are available for removal to
//   grow the maze.
//======================================================================
static void
Add_To_Available(Cell *cell, int *available, int &num_available)
//======================================================================
{
	int i, j;

	// Add edges from cell to the set that are available for removal to
	// grow the maze.

	for ( i = 0 ; i < 4 ; i++ ){
		Cell    *neighbor = cell->edges[i]->Neighbor(cell);

		if ( neighbor && ! neighbor->counter )	{
			int candidate = cell->edges[i]->index;
			for ( j = 0 ; j < num_available ; j++ )
				if ( candidate == available[j] ) {
					printf("Breaking early\n");
					break;
			}
			if ( j == num_available )  {
				available[num_available] = candidate;
				num_available++;
			}
		}
	}

	cell->counter = 1;
}


//**********************************************************************
//
// * Grow a maze by removing candidate edges until all the cells are
//   connected. The edges are not actually removed, they are just made
//   transparent.
//======================================================================
void Maze::
Build_Maze()
//======================================================================
{
	Cell    *to_expand;
	int     index;
	int     *available = new int[num_edges];
	int     num_available = 0;
	int	    num_visited;
	int	    i;

	srand(time(NULL));

	// Choose a random starting cell.
	index = (int)floor((rand() / (float)RAND_MAX) * num_cells);
	to_expand = cells[index];
	Add_To_Available(to_expand, available, num_available);
	num_visited = 1;

	// Join cells up by making edges opaque.
	while ( num_visited < num_cells && num_available > 0 ) {
		int ei;

		index = (int)floor((rand() / (float)RAND_MAX) * num_available);
		to_expand = NULL;

		ei = available[index];

		if ( edges[ei]->neighbors[0] && 
			 !edges[ei]->neighbors[0]->counter )
			to_expand = edges[ei]->neighbors[0];
		else if ( edges[ei]->neighbors[1] && 
			 !edges[ei]->neighbors[1]->counter )
			to_expand = edges[ei]->neighbors[1];

		if ( to_expand ) {
			edges[ei]->opaque = false;
			Add_To_Available(to_expand, available, num_available);
			num_visited++;
		}

		available[index] = available[num_available-1];
		num_available--;
	}

	for ( i = 0 ; i < num_cells ; i++ )
		cells[i]->counter = 0;
}


//**********************************************************************
//
// * Go through all the vertices looking for the minimum and maximum
//   extents of the maze.
//======================================================================
void Maze::
Set_Extents(void)
//======================================================================
{
	int i;

	min_xp = vertices[0]->posn[Vertex::X];
	max_xp = vertices[0]->posn[Vertex::X];
	min_yp = vertices[0]->posn[Vertex::Y];
	max_yp = vertices[0]->posn[Vertex::Y];
	for ( i = 1 ; i < num_vertices ; i++ ) {
		if ( vertices[i]->posn[Vertex::X] > max_xp )
			 max_xp = vertices[i]->posn[Vertex::X];
		if ( vertices[i]->posn[Vertex::X] < min_xp )
			 min_xp = vertices[i]->posn[Vertex::X];
		if ( vertices[i]->posn[Vertex::Y] > max_yp )
			 max_yp = vertices[i]->posn[Vertex::Y];
		if ( vertices[i]->posn[Vertex::Y] < min_yp )
			 min_yp = vertices[i]->posn[Vertex::Y];
    }
}


//**********************************************************************
//
// * Figure out which cell the view is in, using seed_cell as an
//   initial guess. This procedure works by repeatedly checking
//   whether the viewpoint is in the current cell. If it is, we're
//   done. If not, Point_In_Cell returns in new_cell the next cell
//   to test. The new cell is the one on the other side of an edge
//   that the point is "outside" (meaning that it might be inside the
//   new cell).
//======================================================================
void Maze::Find_View_Cell(Cell* seed_cell)
//======================================================================
{
	Cell* new_cell;

	// 
	while (!(seed_cell->Point_In_Cell(viewer_posn[X], viewer_posn[Y], viewer_posn[Z], new_cell))) {
		if (new_cell == 0) {
			// The viewer is outside the top or bottom of the maze.
			throw new MazeException("Maze: View not in maze\n");
		}

		seed_cell = new_cell;
	}

	view_cell = seed_cell;
}


//**********************************************************************
//
// * Move the viewer's position. This method will do collision detection
//   between the viewer's location and the walls of the maze and prevent
//   the viewer from passing through walls.
//======================================================================
void Maze::Move_View_Posn(const float dx, const float dy, const float dz)
//======================================================================
{
	Cell* new_cell;
	float   xs, ys, zs, xe, ye, ze;

	// Move the viewer by the given amount. This does collision testing to
	// prevent walking through walls. It also keeps track of which cells the
	// viewer is in.

	// Set up a line segment from the start to end points of the motion.
	xs = viewer_posn[X];
	ys = viewer_posn[Y];
	zs = viewer_posn[Z];
	xe = xs + dx;
	ye = ys + dy;
	ze = zs + dz;

	// Fix the z to keep it in the maze.
	if (ze > 1.0f - BUFFER);//ze = 1.0f - BUFFER;
	if (ze < BUFFER - 1.0f);//ze = BUFFER - 1.0f;

	// Clip_To_Cell clips the motion segment to the view_cell if the
	// segment intersects an opaque edge. If the segment intersects
	// a transparent edge (through which it can pass), then it clips
	// the motion segment so that it _starts_ at the transparent edge,
	// and it returns the cell the viewer is entering. We keep going
	// until Clip_To_Cell returns NULL, meaning we've done as much of
	// the motion as is possible without passing through walls.
	while ((new_cell = view_cell->Clip_To_Cell(xs, ys, xe, ye, BUFFER)))
		view_cell = new_cell;

	// The viewer is at the end of the motion segment, which may have
	// been clipped.
	viewer_posn[X] = xe;
	viewer_posn[Y] = ye;
	viewer_posn[Z] = ze;
}

//**********************************************************************
//
// * Set the viewer's location 
//======================================================================
void Maze::
Set_View_Posn(float x, float y, float z)
//======================================================================
{
	// First make sure it's in some cell.
	// This assumes that the maze is rectangular.
	if ( x < min_xp + BUFFER )
		x = min_xp + BUFFER;
	if ( x > max_xp - BUFFER )
		x = max_xp - BUFFER;
	if ( y < min_yp + BUFFER )
		y = min_yp + BUFFER;
	if ( y > max_yp - BUFFER )
		y = max_yp - BUFFER;
	if ( z < -1.0f + BUFFER )
		z = -1.0f + BUFFER;
	if ( z > 1.0f - BUFFER )
		z = 1.0f - BUFFER;

	viewer_posn[X] = x;
	viewer_posn[Y] = y;
	viewer_posn[Z] = z;

	// Figure out which cell we're in.
	Find_View_Cell(cells[0]);
}


//**********************************************************************
//
// * Set the angle in which the viewer is looking.
//======================================================================
void Maze::
Set_View_Dir(const float d)
//======================================================================
{
	viewer_dir = d;
}


//**********************************************************************
//
// * Set the horizontal field of view.
//======================================================================
void Maze::
Set_View_FOV(const float f)
//======================================================================
{
	viewer_fov = f;
}


//**********************************************************************
//
// * Draws the map view of the maze. It is passed the minimum and maximum
//   corners of the window in which to draw.
//======================================================================
void Maze::
Draw_Map(int min_x, int min_y, int max_x, int max_y)
//======================================================================
{
	int	    height;
	float   scale_x, scale_y, scale;
	int	    i;

	// Figure out scaling factors and the effective height of the window.
	scale_x = ( max_x - min_x - 10 ) / ( max_xp - min_xp );
	scale_y = ( max_y - min_y - 10 ) / ( max_yp - min_yp );
	scale = scale_x > scale_y ? scale_y : scale_x;
	height = (int)ceil(scale * ( max_yp - min_yp ));

	min_x += 5;
	min_y += 5;

	// Draw all the opaque edges.
	for ( i = 0 ; i < num_edges ; i++ )
		if ( edges[i]->opaque )	{
			float   x1, y1, x2, y2;

			x1 = edges[i]->endpoints[Edge::START]->posn[Vertex::X];
			y1 = edges[i]->endpoints[Edge::START]->posn[Vertex::Y];
			x2 = edges[i]->endpoints[Edge::END]->posn[Vertex::X];
			y2 = edges[i]->endpoints[Edge::END]->posn[Vertex::Y];

			fl_color((unsigned char)floor(edges[i]->color[0] * 255.0),
					 (unsigned char)floor(edges[i]->color[1] * 255.0),
					 (unsigned char)floor(edges[i]->color[2] * 255.0));
			fl_line_style(FL_SOLID);
			fl_line(min_x + (int)floor((x1 - min_xp) * scale),
					  min_y + height - (int)floor((y1 - min_yp) * scale),
					  min_x + (int)floor((x2 - min_xp) * scale),
					  min_y + height - (int)floor((y2 - min_yp) * scale));
		}
}


//**********************************************************************
//
// * Draws the first-person view of the maze. It is passed the focal distance.
//   THIS IS THE FUINCTION YOU SHOULD MODIFY.
//======================================================================
void Maze::
Draw_View(const float focal_dist, Matrix4 projection, Matrix4 modelview)
//======================================================================
{
	frame_num++;

	//###################################################################
	// TODO
	// The rest is up to you!
	//###################################################################
	


	//cout << "current cell: " << view_cell->index << endl;

	vector<vector<Vector4>> edgePoints_in_view = clip_edges();
	
	//visibility + clipping
	//vector<vector<Vector4>> edgePoints_in_view = edges_visibility();

	projection_matrix = projection;
	modelview_matrix = modelview;
	

	//Can't use these functions
	//glClear(GL_DEPTH_BUFFER_BIT);
	//glEnable(GL_DEPTH_TEST);
	for (int i = 0; i < edgePoints_in_view.size(); ++i) {

		float color[3] = {
			edges[i]->color[0],
			edges[i]->color[1],
			edges[i]->color[2]
		};

		if (edges[i]->opaque) {
			Draw_Wall(edgePoints_in_view[i], color);
		}
	}

}

//Draw the wall function
void Maze::Draw_Wall(const std::vector<Vector4> vertices, const float color[3]) {

	//four vertices contstruct a wall
	Vector4 edgeBegin1 = vertices[0];
	Vector4 edgeEnd1 = vertices[1];
	Vector4 edgeEnd2 = vertices[2];
	Vector4 edgeBegin2 = vertices[3];

	////mutiply the vertices by the modelview matrix
	//edgeBegin1 = modelview_matrix * edgeBegin1;
	//edgeEnd1 = modelview_matrix * edgeEnd1;
	//edgeEnd2 = modelview_matrix * edgeEnd2;
	//edgeBegin2 = modelview_matrix * edgeBegin2;

	//cout << edgeBegin1 << endl;
	//cout << edgeEnd1 << endl;
	//cout << endl;
	//cout << endl;

	//mutiply the vertices by the projection matrix
	edgeBegin1 = projection_matrix * edgeBegin1;
	edgeEnd1 = projection_matrix * edgeEnd1;
	edgeEnd2 = projection_matrix * edgeEnd2;
	edgeBegin2 = projection_matrix * edgeBegin2;

	//cout << edgeBegin1 << endl;
	//cout << edgeEnd1 << endl;
	//cout << endl;
	//if (edgeBegin1.w < 0 || edgeEnd1.w < 0) return;


	glBegin(GL_POLYGON);
	glColor3fv(color);

	//x, y, z, w given
	//glVertex4f(edgeBegin1.x, edgeBegin1.y, edgeBegin1.z, edgeBegin1.w);
	//glVertex4f(edgeEnd1.x, edgeEnd1.y, edgeEnd1.z, edgeEnd1.w);
	//glVertex4f(edgeEnd2.x, edgeEnd2.y, edgeEnd2.z, edgeEnd2.w);
	//glVertex4f(edgeBegin2.x, edgeBegin2.y, edgeBegin2.z, edgeBegin2.w);

	edgeBegin1 /= abs(edgeBegin1.w);
	edgeEnd1 /= abs(edgeEnd1.w);
	edgeEnd2 /= abs(edgeEnd2.w);
	edgeBegin2 /= abs(edgeBegin2.w);

	glVertex2f(edgeBegin1.x, edgeBegin1.y);
	glVertex2f(edgeEnd1.x, edgeEnd1.y);
	glVertex2f(edgeEnd2.x, edgeEnd2.y);
	glVertex2f(edgeBegin2.x, edgeBegin2.y);

	glEnd();
}

//**********************************************************************
//
// * Draws the frustum on the map view of the maze. It is passed the
//   minimum and maximum corners of the window in which to draw.
//======================================================================
void Maze::
Draw_Frustum(int min_x, int min_y, int max_x, int max_y)
//======================================================================
{
	int	  height;
	float   scale_x, scale_y, scale;
	float   view_x, view_y;

	// Draws the view frustum in the map. Sets up all the same viewing
	// parameters as draw().
	scale_x = (max_x - min_x - 10) / (max_xp - min_xp);
	scale_y = (max_y - min_y - 10) / (max_yp - min_yp);
	scale = scale_x > scale_y ? scale_y : scale_x;
	height = (int)ceil(scale * (max_yp - min_yp));

	min_x += 5;
	min_y += 5;

	view_x = (viewer_posn[X] - min_xp) * scale;
	view_y = (viewer_posn[Y] - min_yp) * scale;

	fl_line(
		min_x + (int)(view_x + cos(To_Radians(viewer_dir + viewer_fov / 2.0)) * scale),
		min_y + height - (int)(view_y + sin(To_Radians(viewer_dir + viewer_fov / 2.0)) * scale),
		min_x + (int)(view_x),
		min_y + height - (int)(view_y));
	//system("cls");
	//cout << "x0: " << min_x + (int)floor(view_x + cos(To_Radians(viewer_dir + viewer_fov / 2.0)) * scale) << endl;
	//cout << "x1: " << min_y + height - (int)floor(view_y + sin(To_Radians(viewer_dir + viewer_fov / 2.0)) * scale) << endl;
	//cout << "y0: " << min_x + (int)floor(view_x) << endl;
	//cout << "y0: " << min_y + height - (int)floor(view_y) << endl;
	//cout << endl;
	fl_line(
		min_x + (int)(view_x + cos(To_Radians(viewer_dir - viewer_fov / 2.0)) * scale),
		min_y + height - (int)(view_y + sin(To_Radians(viewer_dir - viewer_fov / 2.0)) * scale),
		min_x + (int)(view_x),
		min_y + height - (int)(view_y));
}


//**********************************************************************
//
// * Draws the viewer's cell and its neighbors in the map view of the maze.
//   It is passed the minimum and maximum corners of the window in which
//   to draw.
//======================================================================
void Maze::
Draw_Neighbors(int min_x, int min_y, int max_x, int max_y)
//======================================================================
{
	int	    height;
	float   scale_x, scale_y, scale;
	int	    i, j;

	// Draws the view cell and its neighbors in the map. This works
	// by drawing just the neighbor's edges if there is a neighbor,
	// otherwise drawing the edge. Every edge is shared, so drawing the
	// neighbors' edges also draws the view cell's edges.

	scale_x = ( max_x - min_x - 10 ) / ( max_xp - min_xp );
	scale_y = ( max_y - min_y - 10 ) / ( max_yp - min_yp );
	scale = scale_x > scale_y ? scale_y : scale_x;
	height = (int)ceil(scale * ( max_yp - min_yp ));

	min_x += 5;
	min_y += 5;

	for ( i = 0 ; i < 4 ; i++ )   {
		Cell	*neighbor = view_cell->edges[i]->Neighbor(view_cell);

		if ( neighbor ){
			for ( j = 0 ; j < 4 ; j++ ){
				Edge    *e = neighbor->edges[j];

				if ( e->opaque )	{
					float   x1, y1, x2, y2;

					x1 = e->endpoints[Edge::START]->posn[Vertex::X];
					y1 = e->endpoints[Edge::START]->posn[Vertex::Y];
					x2 = e->endpoints[Edge::END]->posn[Vertex::X];
					y2 = e->endpoints[Edge::END]->posn[Vertex::Y];

					fl_color((unsigned char)floor(e->color[0] * 255.0),
							  (unsigned char)floor(e->color[1] * 255.0),
							  (unsigned char)floor(e->color[2] * 255.0));
					fl_line_style(FL_SOLID);
					fl_line( min_x + (int)floor((x1 - min_xp) * scale),
							 min_y + height - (int)floor((y1 - min_yp) * scale),
							 min_x + (int)floor((x2 - min_xp) * scale),
							 min_y + height - (int)floor((y2 - min_yp) * scale));
				}
			}
		}
		else {
			Edge    *e = view_cell->edges[i];

			if ( e->opaque ){
				float   x1, y1, x2, y2;

				x1 = e->endpoints[Edge::START]->posn[Vertex::X];
				y1 = e->endpoints[Edge::START]->posn[Vertex::Y];
				x2 = e->endpoints[Edge::END]->posn[Vertex::X];
				y2 = e->endpoints[Edge::END]->posn[Vertex::Y];

				fl_color((unsigned char)floor(e->color[0] * 255.0),
							 (unsigned char)floor(e->color[1] * 255.0),
							 (unsigned char)floor(e->color[2] * 255.0));
				fl_line_style(FL_SOLID);
				fl_line(min_x + (int)floor((x1 - min_xp) * scale),
							min_y + height - (int)floor((y1 - min_yp) * scale),
							min_x + (int)floor((x2 - min_xp) * scale),
							min_y + height - (int)floor((y2 - min_yp) * scale));
			 }
		}
	}
}


//**********************************************************************
//
// * Save the maze to a file of the given name.
//======================================================================
bool Maze::
Save(const char *filename)
//======================================================================
{
	FILE    *f = fopen(filename, "w");
	int	    i;

	// Dump everything to a file of the given name. Returns false if it
	// couldn't open the file. True otherwise.

	if ( ! f )  {
		return false;
   }

	fprintf(f, "%d\n", num_vertices);
	for ( i = 0 ; i < num_vertices ; i++ )
		fprintf(f, "%g %g\n", vertices[i]->posn[Vertex::X],
			      vertices[i]->posn[Vertex::Y]);

		fprintf(f, "%d\n", num_edges);
	for ( i = 0 ; i < num_edges ; i++ )
	fprintf(f, "%d %d %d %d %d %g %g %g\n",
				edges[i]->endpoints[Edge::START]->index,
				edges[i]->endpoints[Edge::END]->index,
				edges[i]->neighbors[Edge::LEFT] ?
				edges[i]->neighbors[Edge::LEFT]->index : -1,
				edges[i]->neighbors[Edge::RIGHT] ?
				edges[i]->neighbors[Edge::RIGHT]->index : -1,
				edges[i]->opaque ? 1 : 0,
				edges[i]->color[0], edges[i]->color[1], edges[i]->color[2]);

	fprintf(f, "%d\n", num_cells);
	for ( i = 0 ; i < num_cells ; i++ )
		fprintf(f, "%d %d %d %d\n",
					cells[i]->edges[0] ? cells[i]->edges[0]->index : -1,
					cells[i]->edges[1] ? cells[i]->edges[1]->index : -1,
					cells[i]->edges[2] ? cells[i]->edges[2]->index : -1,
					cells[i]->edges[3] ? cells[i]->edges[3]->index : -1);

	   fprintf(f, "%g %g %g %g %g\n",
					viewer_posn[X], viewer_posn[Y], viewer_posn[Z],
					viewer_dir, viewer_fov);

	fclose(f);

	return true;
}

//generate clip plane with given two fov angles
vector<LineSeg> make_clip_lines(float fov_start, float fov_end) {
	//clip plane is formed by these four lines
	//���@�u (in x right, z up coordinate)
	//line1, near
	LineSeg line1(
		0.01 * tan(Maze::To_Radians(fov_start / 2.0)),
		0.01,
		0.01 * tan(Maze::To_Radians(fov_end / 2.0)),
		0.01
	);
	//line2 left
	LineSeg line2(
		0.01 * tan(Maze::To_Radians(fov_end / 2.0)),
		0.01,
		200 * tan(Maze::To_Radians(fov_end / 2.0)),
		200
	);
	//line3 far
	LineSeg line3(
		200 * tan(Maze::To_Radians(fov_end / 2.0)),
		200,
		200 * tan(Maze::To_Radians(fov_start / 2.0)),
		200
	);
	//line4 right
	LineSeg line4(
		200 * tan(Maze::To_Radians(fov_start / 2.0)),
		200,
		0.01 * tan(Maze::To_Radians(fov_start / 2.0)),
		0.01
	);
	return vector<LineSeg>({ line1, line2, line3, line4});
}

//clip an edge with a clip line
void clip(Edge& edge, const LineSeg clip_line) {
	LineSeg edgeLine(
		edge.endpoints[Edge::START]->posn[Vertex::X],
		edge.endpoints[Edge::START]->posn[Vertex::Y],
		edge.endpoints[Edge::END]->posn[Vertex::X],
		edge.endpoints[Edge::END]->posn[Vertex::Y]
	);
	float edge_x1 = edge.endpoints[Edge::START]->posn[Vertex::X];
	float edge_y1 = edge.endpoints[Edge::START]->posn[Vertex::Y];
	float edge_x2 = edge.endpoints[Edge::END]->posn[Vertex::X];
	float edge_y2 = edge.endpoints[Edge::END]->posn[Vertex::Y];

	float edge_vector[2] = {
		edge.endpoints[Edge::END]->posn[Vertex::X] - edge.endpoints[Edge::START]->posn[Vertex::X],
		edge.endpoints[Edge::END]->posn[Vertex::Y] - edge.endpoints[Edge::START]->posn[Vertex::Y]
	};

	Edge clip_edge(0,
		new Vertex(0,clip_line.start[0], clip_line.start[1]), 
		new Vertex(0, clip_line.end[0], clip_line.end[1]),
		1.0, 1.0, 1.0);

	//��edge����ӻ��A�Pclip_line��cross_param > 1 �άO < 0�A�N���̨S�����I
	float cross_param = edgeLine.Cross_Param(clip_line);
	// Case 1 : edge and clip_line has intersection
	if (cross_param > 0 && cross_param < 1) {
		//intersectino point
		Vertex newPoint(0, edgeLine.start[0] + cross_param*edge_vector[0], edgeLine.start[1] + cross_param * edge_vector[1]);

		//�p�Gedge start �b clip line ���k��(����)�A�h�s�Jnew point & edge start
		//�_�h�s�J new point & edge end
		char point_side = clip_edge.Point_Side(edge_x1, edge_y1);
		if (point_side == Edge::RIGHT) {
			edge.endpoints[Edge::END]->posn[Vertex::X] = newPoint.posn[Vertex::X];
			edge.endpoints[Edge::END]->posn[Vertex::Y] = newPoint.posn[Vertex::Y];
		}
		else {
			edge.endpoints[Edge::START]->posn[Vertex::X] = newPoint.posn[Vertex::X];
			edge.endpoints[Edge::START]->posn[Vertex::Y] = newPoint.posn[Vertex::Y];
		}
	}

	// Case 2: edge and clip_line has no intersection
	if (cross_param < 0 || cross_param > 1) {
		//�p�G���@�I �b clip line ���k��(����)�A�h�����edge����
		//�_�h�Nedge�]��-1
		char point_side = clip_edge.Point_Side(edge_x1, edge_y1);
		if (point_side == Edge::RIGHT) {
			//don't need to modify
		}
		else {
			edge.endpoints[Edge::START]->posn[Vertex::X] = -1;
			edge.endpoints[Edge::START]->posn[Vertex::Y] = -1;
			edge.endpoints[Edge::END]->posn[Vertex::X] = -1;
			edge.endpoints[Edge::END]->posn[Vertex::Y] = -1;
		}
	}
}

vector<vector<Vector4>> Maze::clip_edges() {
	//construct output vertices
	vector<vector<Vector4>> output_vertices;


	//four clip lines: near left far right
	vector<LineSeg> clip_lines = make_clip_lines(viewer_fov, -viewer_fov);

	//LineSeg testLine1(0,1,2,1);
	//LineSeg testLine2(-2,0,-2,3);
	//cout << "crossPraram: " << testLine1.Cross_Param(testLine2) << endl;
	//Edge clip_edge(0,
	//	new Vertex(0, testLine1.start[0], testLine1.start[1]),
	//	new Vertex(0, testLine1.end[0], testLine1.end[1]),
	//	1.0, 1.0, 1.0);
	//cout << "point side: " << (int)clip_edge.Point_Side(1, -1) << endl;


	//iterate all the edges, and clip them to the output_edges
	for (int i = 0; i < (int)num_edges; ++i) {
		//cout << "i: " << i << endl;

		//initialize vertices
		float edge_start[2] = {
			edges[i]->endpoints[Edge::START]->posn[Vertex::X],
			edges[i]->endpoints[Edge::START]->posn[Vertex::Y]
		};
		float edge_end[2] = {
			edges[i]->endpoints[Edge::END]->posn[Vertex::X],
			edges[i]->endpoints[Edge::END]->posn[Vertex::Y]
		};
		Vector4 edgeBegin1(edge_start[Y], 1.0, edge_start[X], 1);
		Vector4 edgeEnd1(edge_end[Y], 1.0, edge_end[X], 1);
		Vector4 edgeEnd2(edge_end[Y], -1.0, edge_end[X], 1);
		Vector4 edgeBegin2(edge_start[Y], -1.0, edge_start[X], 1);

		//mutiply the vertices by the modelview matrix
		edgeBegin1 = modelview_matrix * edgeBegin1;
		edgeEnd1 = modelview_matrix * edgeEnd1;
		edgeEnd2 = modelview_matrix * edgeEnd2;
		edgeBegin2 = modelview_matrix * edgeBegin2;

		//after mutiply the modelview matrix, these points are in this format
		//edgestart(x, y)=>(y, height, x, w) vector4(x, y, z, w), the z is negative the the point is in front of camera

		//cout << i << endl;
		//cout << "Before: " << endl;
		//cout << edgeBegin1 << endl;
		//cout << edgeEnd1 << endl;
		//cout << endl;
		//cout << endl;
		
		//mutiply the z value by -1 to make it x right z up coordinate
		//�i��ݭn�Ѻc�A�ݰ_��vertex�S������
		Edge current_edge(
			i,
			new Vertex(0, edgeBegin1.x, -edgeBegin1.z),
			new Vertex(0, edgeEnd1.x, -edgeEnd1.z),
			edges[i]->color[0],
			edges[i]->color[1],
			edges[i]->color[2]
		);

		//clip(current_edge, clip_lines[0]);
		for (int j = 0; j < clip_lines.size(); ++j) {
			clip(current_edge, clip_lines[j]);
			//outside of the clip plane
			if (current_edge.endpoints[Edge::START]->posn[Vertex::X] == -1 &&
				current_edge.endpoints[Edge::START]->posn[Vertex::Y] == -1 &&
				current_edge.endpoints[Edge::END]->posn[Vertex::Y] == -1 &&
				current_edge.endpoints[Edge::END]->posn[Vertex::Y] == -1) {
				break;
			}
		}

		edgeBegin1.x = current_edge.endpoints[Edge::START]->posn[Vertex::X];
		edgeBegin1.z = -current_edge.endpoints[Edge::START]->posn[Vertex::Y];

		edgeEnd1.x = current_edge.endpoints[Edge::END]->posn[Vertex::X];
		edgeEnd1.z = -current_edge.endpoints[Edge::END]->posn[Vertex::Y];

		edgeEnd2.x = current_edge.endpoints[Edge::END]->posn[Vertex::X];
		edgeEnd2.z = -current_edge.endpoints[Edge::END]->posn[Vertex::Y];

		edgeBegin2.x = current_edge.endpoints[Edge::START]->posn[Vertex::X];
		edgeBegin2.z = -current_edge.endpoints[Edge::START]->posn[Vertex::Y];
		
		//cout << i << endl;
		//cout << "After: " << endl;
		//cout << edgeBegin1 << endl;
		//cout << edgeEnd1 << endl;
		//cout << endl;
		//cout << endl;

		//float x0, x1, y0, y1;
		//x0 = edges[i]->endpoints[Edge::START]->posn[Vertex::X];
		//x1 = edges[i]->endpoints[Edge::END]->posn[Vertex::X];
		//y0 = edges[i]->endpoints[Edge::START]->posn[Vertex::Y];
		//y1 = edges[i]->endpoints[Edge::END]->posn[Vertex::Y];
		//float end_points[2][2] = { {x0,y0},{x1,y1} };



		vector<Vector4> current_vertices = {
			edgeBegin1,
			edgeEnd1,
			edgeEnd2,
			edgeBegin2
		};
		output_vertices.push_back(current_vertices);
	}

	return output_vertices;
}

vector<vector<Vector4>> Maze::edges_visibility() {
	vector<Vector4> _init_vector(4, Vector4(-1,-1,-1,-1));
	vector<vector<Vector4>> result_edges(num_edges, _init_vector);
	set<Cell> used_cell;

	recursive_visibility(*view_cell, viewer_fov, -viewer_fov, used_cell, result_edges);

	return result_edges;
}

void Maze::recursive_visibility(
	const Cell& current_cell, 
	const float fov_start, 
	const float fov_end, 
	std::set<Cell>& used_cell, 
	std::vector<std::vector<Vector4>>& vertices) {
	
	//four clip lines: near left far right
	vector<LineSeg> clip_lines = make_clip_lines(fov_start, fov_end);

	//store edges in current cell and in the camera
	vector<Edge> view_edges;
	for (int i = 0; i < 4; ++i) {
		bool is_visible = true;

		//initialize vertices
		float edge_start[2] = {
			current_cell.edges[i]->endpoints[Edge::START]->posn[Vertex::X],
			current_cell.edges[i]->endpoints[Edge::START]->posn[Vertex::Y]
		};
		float edge_end[2] = {
			current_cell.edges[i]->endpoints[Edge::END]->posn[Vertex::X],
			current_cell.edges[i]->endpoints[Edge::END]->posn[Vertex::Y]
		};
		Vector4 edgeBegin1(edge_start[Y], 1.0, edge_start[X], 1);
		Vector4 edgeEnd1(edge_end[Y], 1.0, edge_end[X], 1);
		Vector4 edgeEnd2(edge_end[Y], -1.0, edge_end[X], 1);
		Vector4 edgeBegin2(edge_start[Y], -1.0, edge_start[X], 1);

		//mutiply the vertices by the modelview matrix
		edgeBegin1 = modelview_matrix * edgeBegin1;
		edgeEnd1 = modelview_matrix * edgeEnd1;
		edgeEnd2 = modelview_matrix * edgeEnd2;
		edgeBegin2 = modelview_matrix * edgeBegin2;

		Edge current_edge(
			i,
			new Vertex(0, edgeBegin1.x, -edgeBegin1.z),
			new Vertex(0, edgeEnd1.x, -edgeEnd1.z),
			edges[i]->color[0],
			edges[i]->color[1],
			edges[i]->color[2]
		);

		//clip(current_edge, clip_lines[0]);
		for (int j = 0; j < clip_lines.size(); ++j) {
			clip(current_edge, clip_lines[j]);
			//outside of the clip plane
			if (current_edge.endpoints[Edge::START]->posn[Vertex::X] == -1 &&
				current_edge.endpoints[Edge::START]->posn[Vertex::Y] == -1 &&
				current_edge.endpoints[Edge::END]->posn[Vertex::Y] == -1 &&
				current_edge.endpoints[Edge::END]->posn[Vertex::Y] == -1) {
				is_visible = false;
				break;
			}
		}

		if (is_visible) {
			edgeBegin1.x = current_edge.endpoints[Edge::START]->posn[Vertex::X];
			edgeBegin1.z = -current_edge.endpoints[Edge::START]->posn[Vertex::Y];

			edgeEnd1.x = current_edge.endpoints[Edge::END]->posn[Vertex::X];
			edgeEnd1.z = -current_edge.endpoints[Edge::END]->posn[Vertex::Y];

			edgeEnd2.x = current_edge.endpoints[Edge::END]->posn[Vertex::X];
			edgeEnd2.z = -current_edge.endpoints[Edge::END]->posn[Vertex::Y];

			edgeBegin2.x = current_edge.endpoints[Edge::START]->posn[Vertex::X];
			edgeBegin2.z = -current_edge.endpoints[Edge::START]->posn[Vertex::Y];
		}




	}
}