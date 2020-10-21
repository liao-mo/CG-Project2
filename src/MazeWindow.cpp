/************************************************************************
     File:        MazeWindow.cpp

     Author:     
                  Stephen Chenney, schenney@cs.wisc.edu
     Modifier
                  Yu-Chi Lai, yu-chi@cs.wisc.edu

     Comment:    
						(c) 2001-2002 Stephen Chenney, University of Wisconsin at Madison

						Class header file for the MazeWindow class. The MazeWindow is
						the window in which the viewer's view of the maze is displayed.
		

     Platform:    Visio Studio.Net 2003 (converted to 2005)

*************************************************************************/

#include "MazeWindow.h"
#include <Fl/math.h>
#include <Fl/gl.h>
#include <GL/glu.h>
#include <stdio.h>
#include "Matrices.h"
#include <iostream>
#include <iomanip>

using namespace std;


//*************************************************************************
//
// * Constructor
//=========================================================================
MazeWindow::
MazeWindow(int x, int y, int width, int height, const char *label,Maze *m)
	: Fl_Gl_Window(x, y, width, height, label)
//=========================================================================
{
	maze = m;

	// The mouse button isn't down and there is no key pressed.
	down = false;
	z_key = 0;
}


//*************************************************************************
//
// * Set the maze. Also causes a redraw.
//=========================================================================
void MazeWindow::
Set_Maze(Maze *m)
//=========================================================================
{
	// Change the maze
	maze = m;

	// Force a redraw
	redraw();
}


//*************************************************************************
//
// * draw() method invoked whenever the view changes or the window
//   otherwise needs to be redrawn.
//=========================================================================
void MazeWindow::
draw(void)
//=========================================================================
{
	float   focal_length;

	if (!valid()) {
		// The OpenGL context may have been changed
		// Set up the viewport to fill the window.
		glViewport(0, 0, w(), h());

		// We are using orthogonal viewing for 2D. This puts 0,0 in the
		// middle of the screen, and makes the image size in view space
		// the same size as the window.
		gluOrtho2D(-w() * 0.5, w() * 0.5, -h() * 0.5, h() * 0.5);

		// Sets the clear color to black.
		glClearColor(0.0, 0.0, 0.0, 1.0);
	}

	// Clear the screen.
	glClear(GL_COLOR_BUFFER_BIT);

	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
	
	glBegin(GL_QUADS);
	// Draw the "floor". It is an infinite plane perpendicular to
	// vertical, so we know it projects to cover the entire bottom
	// half of the screen. Walls of the maze will be drawn over the top
	// of it.
	glColor3f(0.2f, 0.2f, 0.2f);
	glVertex2f(-w() * 0.5f, -h() * 0.5f);
	glVertex2f(w() * 0.5f, -h() * 0.5f);
	glVertex2f(w() * 0.5f, 0.0);
	glVertex2f(-w() * 0.5f, 0.0);

	// Draw the ceiling. It will project to the entire top half
	// of the window.
	glColor3f(0.4f, 0.4f, 0.4f);
	glVertex2f(w() * 0.5f, h() * 0.5f);
	glVertex2f(-w() * 0.5f, h() * 0.5f);
	glVertex2f(-w() * 0.5f, 0.0);
	glVertex2f(w() * 0.5f, 0.0);
	glEnd();


	if (maze) {
		// Set the focal length. We can do this because we know the
		// field of view and the size of the image in view space. Note
		// the static member function of the Maze class for converting
		// radians to degrees. There is also one defined for going backwards.
		focal_length = w() / (float)(2.0 * tan(Maze::To_Radians(maze->viewer_fov) * 0.5));

		// Draw the 3D view of the maze (the visible walls.) You write this.
		// Note that all the information that is required to do the
		// transformations and projection is contained in the Maze class,
		// plus the focal length.

		glClear(GL_DEPTH_BUFFER_BIT);

		////////////////////////////
		//build my modelview matrix
		////////////////////////////
		glMatrixMode(GL_PROJECTION);
		glLoadIdentity();

		//can't use this function
		float aspect = (float)(w() / h());
		//gluPerspective(maze->viewer_fov, aspect, 0.01, 200);

		////////////////////////////
		//build my modelview matrix
		////////////////////////////
		Matrix4 projection_matrix = construct_perspective_projection_matrix(maze->viewer_fov, (float)(w() / h()), 0.01, 200);

		//debug
		//display my projection matrix
		//cout << "projection matrix: " << endl;
		//cout << projection_matrix << endl;

		//apply this projection matrix
		//glMultMatrixf(projection_matrix.get());

		//Debug
		//get projection_matrix
		//double _projection_matrix[16];
		//glGetDoublev(GL_PROJECTION_MATRIX, _projection_matrix);
		//cout << "projection matrix 10x10-100: " << endl;
		//for (int i = 0; i < 16; ++i) {
		//	if (i % 4 == 0) cout << endl;
		//	cout << setw(10) << _projection_matrix[i] << " ";
		//}
		//cout << "\n\n";



		////////////////////////////
		//set model view matrix
		////////////////////////////
		glMatrixMode(GL_MODELVIEW);
		glLoadIdentity();

		float viewer_pos[3] = { maze->viewer_posn[Maze::Y], 0.0f, maze->viewer_posn[Maze::X] };
		////////////////////////////
		//can't use this function
		////////////////////////////
		//gluLookAt(
		//	viewer_pos[Maze::X],
		//	viewer_pos[Maze::Y],
		//	viewer_pos[Maze::Z],
		//	viewer_pos[Maze::X] + sin(Maze::To_Radians(maze->viewer_dir)),
		//	viewer_pos[Maze::Y],
		//	viewer_pos[Maze::Z] + cos(Maze::To_Radians(maze->viewer_dir)),
		//	0.0,
		//	1.0,
		//	0.0);

		////////////////////////////
		//build my modelview matrix
		////////////////////////////
		//construct eye position vector
		Vector3 eye_pos;
		eye_pos.set(viewer_pos[Maze::X], viewer_pos[Maze::Y], viewer_pos[Maze::Z]);

		//construct target vector
		Vector3 targetVector;
		targetVector.set(viewer_pos[Maze::X] + sin(Maze::To_Radians(maze->viewer_dir)), viewer_pos[Maze::Y], viewer_pos[Maze::Z] + cos(Maze::To_Radians(maze->viewer_dir)));

		//construct up direction vector
		Vector3 upDirection;
		upDirection.set(0.0, 0.1, 0.0);

		//initialize 4x4 modelview matrix
		Matrix4 modelview_matrix = construct_modelview_matrix(eye_pos, targetVector, upDirection);

		//apply this modelview matrix
		//glMultMatrixf(modelview_matrix.get());

		//cout << "my modelview matrix: " << endl;
		//cout << modelview_matrix << endl;
		

		//Debug
		//get projection_matrix
		//double _modelview_matrix[16];
		//glGetDoublev(GL_MODELVIEW_MATRIX, _modelview_matrix);
		//cout << "glu modelview matrix 10x10-100: " << endl;
		//for (int i = 0; i < 16; ++i) {
		//	if (i % 4 == 0) cout << endl;
		//	cout << setw(10) << _modelview_matrix[i] << " ";
		//}
		//cout << "\n\n";

		maze->Draw_View(focal_length, projection_matrix, modelview_matrix);
	}
}


//*************************************************************************
//
// *
//=========================================================================
bool MazeWindow::
Drag(float dt)
//=========================================================================
{
	float   x_move, y_move, z_move;

	if ( down ) {
		int	dx = x_down - x_last;
		int   dy = y_down - y_last;
		float dist;

		// Set the viewing direction based on horizontal mouse motion.
		maze->Set_View_Dir(d_down + 360.0f * dx / (float)w());

		// Set the viewer's linear motion based on a speed (derived from
		// vertical mouse motion), the elapsed time and the viewing direction.
		dist = 10.0f * dt * dy / (float)h();
		x_move = dist * (float)cos(Maze::To_Radians(maze->viewer_dir));
		y_move = dist * (float)sin(Maze::To_Radians(maze->viewer_dir));
	}
	else {
		x_move = 0.0;
		y_move = 0.0;
	}

	// Update the z posn
	z_move = z_key * 0.1f;
	z_key = 0;

	// Tell the maze how much the view has moved. It may restrict the motion
	// if it tries to go through walls.
	maze->Move_View_Posn(x_move*10, y_move*10, z_move*10);

	return true;
}


//*************************************************************************
//
// *
//=========================================================================
bool MazeWindow::
Update(float dt)
//=========================================================================
{
	// Update the view

	if ( down || z_key ) // Only do anything if the mouse button is down.
		return Drag(dt);

	// Nothing changed, so no need for a redraw.
	return false;
}


//*************************************************************************
//
// *
//=========================================================================
int MazeWindow::
handle(int event)
//=========================================================================
{
	if (!maze)
		return Fl_Gl_Window::handle(event);

	// Event handling routine.
	switch ( event ) {
		case FL_PUSH:
			down = true;
			x_last = x_down = Fl::event_x();
			y_last = y_down = Fl::event_y();
			d_down = maze->viewer_dir;
			return 1;
		case FL_DRAG:
			x_last = Fl::event_x();
			y_last = Fl::event_y();
			return 1;
			case FL_RELEASE:
			down = false;
			return 1;
		case FL_KEYBOARD:
			if ( Fl::event_key() == FL_Up )	{
				z_key = 1;
				return 1;
			}
			if ( Fl::event_key() == FL_Down ){
				z_key = -1;
				return 1;
			}
			if (Fl::event_key() == FL_Left) {
				//z_key = 1;
				return 1;
			}
			if (Fl::event_key() == FL_Right) {
				//z_key = -1;
				return 1;
			}
			return Fl_Gl_Window::handle(event);
		case FL_FOCUS:
		case FL_UNFOCUS:
			return 1;
	}

	// Pass any other event types on the superclass.
	return Fl_Gl_Window::handle(event);
}

Matrix4& MazeWindow::construct_perspective_projection_matrix(float fov_angle, float aspect, float z_near, float z_far) {
	//find the values that can describe the pyramid
	float x_right, x_left, y_top, y_bottom;
	y_top = z_near * tanf(Maze::To_Radians(fov_angle) / 2.0);
	y_bottom = -y_top;
	x_right = y_top * aspect;
	x_left = -x_right;

	float x_width = x_right - x_left;
	float y_height = y_top - y_bottom;

	//construct 4x4 projection matrix
	Matrix4 projection_matrix;
	projection_matrix[0] = 2.0 * z_near / x_width;
	projection_matrix[1] = 0;
	projection_matrix[2] = 0;
	projection_matrix[3] = 0;
	projection_matrix[4] = 0;
	projection_matrix[5] = 2.0 * z_near / y_height;
	projection_matrix[6] = 0;
	projection_matrix[7] = 0;
	projection_matrix[8] = (x_right + x_left) / x_width;
	projection_matrix[9] = (y_top + y_bottom) / y_height;
	projection_matrix[10] = -(z_far + z_near) / (z_far - z_near);
	projection_matrix[11] = -1.0;
	projection_matrix[12] = 0;
	projection_matrix[13] = 0;
	projection_matrix[14] = -2.0 * z_far * z_near / (z_far - z_near);
	projection_matrix[15] = 0;
	
	return projection_matrix;
}

Matrix4& MazeWindow::construct_modelview_matrix(Vector3 eye_pos, Vector3 target_vec, Vector3 up_dir) {
	Matrix4 modelview_matrix;
	modelview_matrix.identity();

	//forwardVector is the vector from eye to target
	Vector3 forwardVector = target_vec - eye_pos;
	forwardVector.normalize(); // make it unit length

	//compute the left vector
	Vector3 leftVector = up_dir.cross(forwardVector); // cross product
	leftVector.normalize(); // make it unit length

	//recompute the orthonormal up vector by forward and left vector
	Vector3 upVector = forwardVector.cross(leftVector); // cross product

	// set rotation part of this matrix
	modelview_matrix[0] = leftVector.x;
	modelview_matrix[4] = leftVector.y;
	modelview_matrix[8] = leftVector.z;
	modelview_matrix[1] = upVector.x;
	modelview_matrix[5] = upVector.y;
	modelview_matrix[9] = upVector.z;
	modelview_matrix[2] = forwardVector.x;
	modelview_matrix[6] = forwardVector.y;
	modelview_matrix[10] = forwardVector.z;

	// set translation part of this matrix
	modelview_matrix[12] = -leftVector.x * eye_pos.x - leftVector.y * eye_pos.y - leftVector.z * eye_pos.z;
	modelview_matrix[13] = -upVector.x * eye_pos.x - upVector.y * eye_pos.y - upVector.z * eye_pos.z;
	modelview_matrix[14] = -forwardVector.x * eye_pos.x - forwardVector.y * eye_pos.y - forwardVector.z * eye_pos.z;

	//mutiply -1 to the matrix expect homogeneous coordinate
	modelview_matrix = -1 * modelview_matrix;
	modelview_matrix[15] = 1;

	return modelview_matrix;
}