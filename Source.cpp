#define OLC_PGE_APPLICATION
#include "olcPixelGameEngine.h"
#include <vector>
#include <iostream>
#define _USE_MATH_DEFINES // for C++
#include <math.h>
using namespace std;
class Donut : public olc::PixelGameEngine
{
private:
	bool play = true; // used to check if paused
	float delphi = 0.01, deltheta = 0.01; // theta and phi angle density
	float alpha1 = 0, alpha2 = 0; // angles of rotation around OX and OY
	float R = 8, r = 5; // radiuses of torus
	float K1 = 300, K2 = 70; // distance from viewer to the screen, distance from torus to screen
	float lx = -1.0, ly = 0, lz = -1; // direction of light
	float sqrt2 = sqrt(2.0);
	int **output; //keeps output pixel saturation (of a red color)
	float **zbuffer; // keeps 1/z of a point projected onto pixel
	int mx, my; // mouse x and y
	float sensitivity = 0.01; // sensitivity of a mouse
public:
	Donut()
	{
		// Name your application
		sAppName = "Donut";
	}
public:
	bool OnUserCreate() override
	{
		// Called once at the start, so create things here
		output = new int*[ScreenWidth()];
		for (int i = 0; i < ScreenWidth(); i++)output[i] = new int[ScreenHeight()];

		zbuffer = new float* [ScreenWidth()];
		for (int i = 0; i < ScreenWidth(); i++)zbuffer[i] = new float[ScreenHeight()];

		for (int i = 0; i < ScreenWidth(); i++)
		{
			for (int j = 0; j < ScreenHeight(); j++)
			{
				output[i][j] = 0;
				zbuffer[i][j] = 0;
			}
		}

		return true;
	}

	bool OnUserUpdate(float fElapsedTime) override
	{
		if (GetMouse(0).bPressed)play = !play; // if left clicked 
		if (!play)		 						//pause animation
		{
			if (GetMouse(1).bPressed) {mx = GetMouseX(); my = GetMouseY();}

			int new_mx, new_my;

			if (GetMouse(1).bHeld)
			{
				new_mx = GetMouseX(); 
				new_my = GetMouseY();

				alpha1 -= sensitivity * (new_mx - mx); // rotate around OX
				alpha2 -= sensitivity * (new_my - my); // rotate around OY

				mx = new_mx; my = new_my;
			}
			else return true;
		}
		else // if animation is not paused
		{
			alpha1 += 1 * fElapsedTime; // auto rotate around OX
			alpha2 += 1 * fElapsedTime; // auto rotate around OY
		}
		Clear(olc::BLACK);

		float cosphi;
		float sinphi;
		float costheta;
		float sintheta;
		float cosalpha1;
		float sinalpha1;
		float cosalpha2;
		float sinalpha2;
		float x, y, z; 
		float tx, ty, tz; // inbetween values 
		float nx, ny, nz; // coordinates of a normal vector
		float tnx, tny, tnz;// inbetween values 
		int xp, yp; // keeps projection coordinates of a point to a screen 
		float ooz; // 1/z
		float res;

		for (float theta = 0; theta < 2 * M_PI; theta += deltheta)
		{
			for (float phi = 0; phi < 2 * M_PI; phi += delphi)
			{
				cosphi = cos(phi);
				sinphi = sin(phi);
				costheta = cos(theta);
				sintheta = sin(theta);
				cosalpha1 = cos(alpha1);
				sinalpha1 = sin(alpha1);
				cosalpha2 = cos(alpha2);
				sinalpha2 = sin(alpha2);

				//calculate x, y, z of point
				x = (R + r * costheta) * cosphi;
				y = (R + r * costheta) * sinphi;
				z = r * sintheta;
				//rotate around OY
				tx = x * cosalpha1 - z * sinalpha1;
				ty = y;
				tz = x * sinalpha1 + z * cosalpha1;
				//rotate around OX
				x = tx;
				y = ty * cosalpha2 - tz * sinalpha2;
				z = ty * sinalpha2 + tz * cosalpha2 + K2;

				//calculate normal vector coordinates
				nx = costheta * cosphi;
				ny = costheta * sinphi;
				nz = sintheta;
				//rotate around OY
				tnx = nx * cosalpha1 - nz * sinalpha1;
				tny = ny;
				tnz = nx * sinalpha1 + nz * cosalpha1;
				//rotate around OX
				nx = tnx;
				ny = tny * cosalpha2 - tnz * sinalpha2;
				nz = tny * sinalpha2 + tnz * cosalpha2;

				ooz = 1 / z;
				//calculate projection on screen coordinates
				xp = (int)(ScreenWidth() / 2 + K1 * ooz * x);
				yp = (int)(ScreenHeight() / 2 - K1 * ooz * y);

				if (xp < 0 || xp >= ScreenWidth())continue;
				if (yp < 0 || yp >= ScreenHeight())continue;
				
				if (zbuffer[xp][yp] < ooz) // if this point is closer to the screen
				{
					zbuffer[xp][yp] = ooz;
					res = nx * lx + ny * ly + nz * lz; // stores cos * sqrt(2) between light direction and normal vector
					res = (res / sqrt2 + 1.0) / 2.0; // make 0 <= res <= 1
					res = pow(100, res) / 100; // do that to make light more drammatic and noticable 
					output[xp][yp] = res * 255;
				}
			}
		}
		//Draw pixels
		for (int i = 0; i < ScreenWidth(); i++)
		{
			for (int j = 0; j < ScreenHeight(); j++)
			{
				Draw(i, j, olc::Pixel(output[i][j], 0, 0));
				zbuffer[i][j] = 0;
				output[i][j] = 0;
			}
		}
		return true;
	}
};

int main()
{
	Donut demo;
	if (demo.Construct(256 * 1, 256 * 1, 3, 3))
		demo.Start();
	return 0;
}