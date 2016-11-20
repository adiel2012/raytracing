/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package reaytracing;

/**
 *
 * @author adiel
 */
public class Sphere
{
    public Sphere(double x, double y, double z, double r, double clr,
                  double clg, double clb)
    {
        cx = x;
        cy = y;
        cz = z;
        radius = r;
        clR = clr;
        clG = clg;
        clB = clb;
    }
    public static double GetCoord(double i1, double i2, double w1, double w2,
        double p)
    {
        return ((p - i1) / (i2 - i1)) * (w2 - w1) + w1;
    }
    public static double modv(double vx, double vy, double vz)
    {
        return Math.sqrt(vx * vx + vy * vy + vz * vz);
    }
    void Move(double vx, double vy, double vz)
    {
        cx += vx;
        cy += vy;
        cz += vz;
    }
    void MoveTo(double vx, double vy, double vz)
    {
        cx = vx;
        cy = vy;
        cz = vz;
    }
    void RotX(double angle)
    {
        double y = cy * Math.cos(angle) - cz * Math.sin(angle);
        double z = cy * Math.sin(angle) + cz * Math.cos(angle);
        cy = y;
        cz = z;
    }
    void RotY(double angle)
    {
        double x = cx * Math.cos(angle) - cz * Math.sin(angle);
        double z = cx * Math.sin(angle) + cz * Math.cos(angle);
        cx = x;
        cz = z;
    }
    public static double GetSphereIntersec(double cx, double cy, double cz,
                         double radius, double px, double py, double pz,
                         double vx, double vy, double vz)
    {
        // x-xo 2 + y-yo 2 + z-zo 2 = r 2
        // x,y,z = p+tv
        // At2 + Bt + C = 0
        double A = (vx * vx + vy * vy + vz * vz);
        double B = 2.0 * (px * vx + py * vy + pz * vz - vx * cx - vy *
                   cy - vz * cz);
        double C = px * px - 2 * px * cx + cx * cx + py * py - 2 * py *
                   cy + cy * cy + pz * pz - 2 * pz * cz + cz * cz -
                   radius * radius;
        double D = B * B - 4 * A * C;
        double t = -1.0;
        if (D >= 0)
        {
            double t1 = (-B - Math.sqrt(D)) / (2.0 * A);
            double t2 = (-B + Math.sqrt(D)) / (2.0 * A);
            if (t1 > t2) t = t1; else t = t2;
        }
        return t;
    }
    public static double GetCosAngleV1V2(double v1x, double v1y, double v1z,
                                         double v2x, double v2y, double v2z)
    {
        /* incident angle
         intersection pt (i)
        double ix, iy, iz;
        ix = px+t*vx;
        iy = py+t*vy;
        iz = pz+t*vz;
        normal at i
        double nx, ny, nz;
        nx = ix - cx;
        ny = iy - cy;
        nz = iz - cz;

        cos(t) = (v.w) / (|v|.|w|)
        */
        return (v1x * v2x + v1y * v2y + v1z * v2z) / (modv(v1x, v1y, v1z) *
               modv(v2x, v2y, v2z));
    }
    public double cx, cy, cz, radius, clR, clG, clB;
}