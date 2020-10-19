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
	if (ze > 1.0f - BUFFER) ze = 1.0f - BUFFER;
	if (ze < BUFFER - 1.0f) ze = BUFFER - 1.0f;

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
	
	vector<Edge> edges_in_view = clip_edges();


	projection_matrix = projection;
	modelview_matrix = modelview;
	

	//Can't use these functions
	//glClear(GL_DEPTH_BUFFER_BIT);
	//glEnable(GL_DEPTH_TEST);
	for (int i = 0; i < (int)num_edges; ++i) {
		float edge_start[2] = {
			edges[i]->endpoints[Edge::START]->posn[Vertex::X],
			edges[i]->endpoints[Edge::START]->posn[Vertex::Y]
		};
		float edge_end[2] = {
			edges[i]->endpoints[Edge::END]->posn[Vertex::X],
			edges[i]->endpoints[Edge::END]->posn[Vertex::Y]
		};
		float color[3] = {
			edges[i]->color[0],
			edges[i]->color[1],
			edges[i]->color[2]
		};

		if (edges[i]->opaque) {
			Draw_Wall(edge_start, edge_end, color);
		}
	}

}

//Draw the wall function
void Maze::Draw_Wall(const float start[2], const float end[2], const float color[3]) {
	//float edge0[3] = { start[Y], 0.0f, start[X] };
	//float edge1[3] = { end[Y], 0.0f, end[X] };

	//four vertices contstruct a wall
	Vector4 edgeBegin1(start[Y], 1.0, start[X], 1);
	Vector4 edgeEnd1(end[Y], 1.0, end[X], 1);
	Vector4 edgeEnd2(end[Y], -1.0, end[X], 1);
	Vector4 edgeBegin2(start[Y], -1.0, start[X], 1);

	//mutiply the vertices by the modelview matrix
	edgeBegin1 = modelview_matrix * edgeBegin1;
	edgeEnd1 = modelview_matrix * edgeEnd1;
	edgeEnd2 = modelview_matrix * edgeEnd2;
	edgeBegin2 = modelview_matrix * edgeBegin2;

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
	//glVertex3f(edge0[X], 1.0f, edge0[Z]);
	//glVertex3f(edge1[X], 1.0f, edge1[Z]);
	//glVertex3f(edge1[X], -1.0f, edge1[Z]);
	//glVertex3f(edge0[X], -1.0f, edge0[Z]);

	//x, y, z, w given
	//glVertex4f(edgeBegin1.x, edgeBegin1.y, edgeBegin1.z, edgeBegin1.w);
	//glVertex4f(edgeEnd1.x, edgeEnd1.y, edgeEnd1.z, edgeEnd1.w);
	//glVertex4f(edgeEnd2.x, edgeEnd2.y, edgeEnd2.z, edgeEnd2.w);
	//glVertex4f(edgeBegin2.x, edgeBegin2.y, edgeBegin2.z, edgeBegin2.w);

	edgeBegin1 /= edgeBegin1.w;
	edgeEnd1 /= edgeEnd1.w;
	edgeEnd2 /= edgeEnd2.w;
	edgeBegin2 /= edgeBegin2.w;

	glVertex3f(edgeBegin1.x, edgeBegin1.y, edgeBegin1.z);
	glVertex3f(edgeEnd1.x, edgeEnd1.y, edgeEnd1.z);
	glVertex3f(edgeEnd2.x, edgeEnd2.y, edgeEnd2.z);
	glVertex3f(edgeBegin2.x, edgeBegin2.y, edgeBegin2.z);

	


	//glVertex2f(edgeBegin1.x, edgeBegin1.y);
	//glVertex2f(edgeEnd1.x, edgeEnd1.y);
	//glVertex2f(edgeEnd2.x, edgeEnd2.y);
	//glVertex2f(edgeBegin2.x, edgeBegin2.y);

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


//clipping functions
float x_intersect(int x1, int y1, int x2, int y2, int x3, int y3, int x4, int y4) {
	float num = (x1 * y2 - y1 * x2) * (x3 - x4) -
		(x1 - x2) * (x3 * y4 - y3 * x4);
	float den = (x1 - x2) * (y3 - y4) - (y1 - y2) * (x3 - x4);
	return num / den;
}
float y_intersect(int x1, int y1, int x2, int y2, int x3, int y3, int x4, int y4) {
	float num = (x1 * y2 - y1 * x2) * (y3 - y4) -
		(y1 - y2) * (x3 * y4 - y3 * x4);
	float den = (x1 - x2) * (y3 - y4) - (y1 - y2) * (x3 - x4);
	return num / den;
}

void clip(float poly_points[][2], int& poly_size, float x1, float y1, float x2, float y2) {
	float new_points[2][2];
	int new_poly_size = 0;

	//(ix,iy),(kx,ky) are the coordinate values of the points 
	//for (int i = 0; i < poly_size; i++)
	//{
		// i and k form a line in polygon 
	int i = 0;
	int k = (i + 1);
	float ix = poly_points[i][0], iy = poly_points[i][1];
	float kx = poly_points[k][0], ky = poly_points[k][1];

	// Calculating position of first point 
	// w.r.t. clipper line 
	float i_pos = (x2 - x1) * (iy - y1) - (y2 - y1) * (ix - x1);

	// Calculating position of second point 
	// w.r.t. clipper line 
	float k_pos = (x2 - x1) * (ky - y1) - (y2 - y1) * (kx - x1);

	// Case 1 : When both points are inside 
	if (i_pos < 0 && k_pos < 0)
	{
		cout << "case1" << endl;
		//Only second point is added 
		new_points[new_poly_size][0] = kx;
		new_points[new_poly_size][1] = ky;
		new_poly_size++;
	}

	// Case 2: When only first point is outside 
	else if (i_pos >= 0 && k_pos < 0)
	{
		cout << "case2" << endl;
		// Point of intersection with edge 
		// and the second point is added 
		new_points[new_poly_size][0] = x_intersect(x1, y1, x2, y2, ix, iy, kx, ky);
		new_points[new_poly_size][1] = y_intersect(x1, y1, x2, y2, ix, iy, kx, ky);
		new_poly_size++;

		new_points[new_poly_size][0] = kx;
		new_points[new_poly_size][1] = ky;
		new_poly_size++;
	}

	// Case 3: When only second point is outside 
	else if (i_pos < 0 && k_pos >= 0)
	{
		cout << "case3" << endl;
		//Only point of intersection with edge is added 
		new_points[new_poly_size][0] = x_intersect(x1, y1, x2, y2, ix, iy, kx, ky);
		new_points[new_poly_size][1] = y_intersect(x1, y1, x2, y2, ix, iy, kx, ky);
		new_poly_size++;
	}

	//// Case 4: When both points are outside 
	else
	{
		cout << "case4" << endl;
		//No points are added 
	}
//}

// Copying new points into original array 
// and changing the no. of vertices 
	poly_size = new_poly_size;
	for (int i = 0; i < poly_size; i++)
	{
		poly_points[i][0] = new_points[i][0];
		poly_points[i][1] = new_points[i][1];
	}
}

vector<Edge> Maze::clip_edges() {
	vector<Edge> output_edges;

	//clip plane is formed by these two lines
	//µøÀ@½u (in x-right, y-up coordinate)
	//line1:composed by 
	//	x1:(view_x)
	//	x2:(view_x + cos(To_Radians(viewer_dir + viewer_fov / 2.0)))
	//	y1:(view_y))
	//	y2:(view_y + sin(To_Radians(viewer_dir + viewer_fov / 2.0)))
	//line2:composed by
	//	x1:(view_x)
	//	x2:(view_x + cos(To_Radians(viewer_dir - viewer_fov / 2.0)))
	//	y1:(view_y))
	//	y2:(view_y + sin(To_Radians(viewer_dir - viewer_fov / 2.0)))

	//line1
	vector<float> line1(4, 0);//(x1,y1,x2,y2)
	line1[0] = viewer_posn[X];
	line1[1] = viewer_posn[Y];
	line1[2] = viewer_posn[X] + cos(To_Radians(viewer_dir + viewer_fov / 2));
	line1[3] = viewer_posn[Y] + sin(To_Radians(viewer_dir + viewer_fov / 2));
	//line2
	vector<float> line2(4, 0);//(x1,y1,x2,y2)
	line2[0] = viewer_posn[X] + cos(To_Radians(viewer_dir - viewer_fov / 2));
	line2[1] = viewer_posn[Y] + sin(To_Radians(viewer_dir - viewer_fov / 2));
	line2[2] = viewer_posn[X];
	line2[3] = viewer_posn[Y];



	//iterate all the edges, and clip them to the output_edges
	for (int i = 0; i < (int)num_edges; ++i) {
		cout << "i: " << i << endl;
		float x0, x1, y0, y1;
		x0 = edges[i]->endpoints[Edge::START]->posn[Vertex::X];
		x1 = edges[i]->endpoints[Edge::END]->posn[Vertex::X];
		y0 = edges[i]->endpoints[Edge::START]->posn[Vertex::Y];
		y1 = edges[i]->endpoints[Edge::END]->posn[Vertex::Y];
		float end_points[2][2] = { {x0,y0},{x1,y1} };
		//clip line1
		int point_size = 2;
		clip(end_points, point_size, line1[0], line1[1], line1[2], line1[3]);
		//clip line2
		if (point_size == 2) {
			//clip(end_points, point_size, line2[0], line2[1], line2[2], line2[3]);
		}

		cout << endl;







		float color[3] = {
			edges[i]->color[0],
			edges[i]->color[1],
			edges[i]->color[2]
		};

		Vertex* temp_v1 = new Vertex (i * 2 + 0, end_points[0][X], end_points[0][Y]);
		Vertex* temp_v2 = new Vertex (i * 2 + 1, end_points[1][X], end_points[1][Y]);
		Edge temp_edge(i, temp_v1, temp_v2, color[0], color[1], color[2]);
		output_edges.push_back(temp_edge);
	}

	return output_edges;
}