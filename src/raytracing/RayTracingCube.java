/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package raytracing;

import java.awt.Color;
import java.awt.Graphics;
import java.awt.Rectangle;
import java.awt.image.BufferedImage;
import java.util.ArrayList;
import java.util.logging.Level;
import java.util.logging.Logger;
import javax.swing.JFrame;

/**
 *
 * @author adiel
 */
public class RayTracingCube {

    public static int W = 400, H = 400;
    public static Lienzo form = new Lienzo();
    // viewer position
    public static double px = (double) 1000,
            py = (double) 1000,
            pz = (double) 1000;
    /**
     * @param args the command line arguments
     */
    //   http://www.codeproject.com/Articles/19732/Simple-Ray-Tracing-in-C

    public static BufferedImage bufferedImage = null;

    public static void main(String[] args) {
        bufferedImage = generarImagen();
        show();
        while (true) {
            bufferedImage = generarImagen();
            form.repaint();
            px -= 5;
            py -= 5;
            pz -= 5;

//            try {
//                Thread.sleep(2000);
//            } catch (InterruptedException ex) {
//                Logger.getLogger(RayTracingCube.class.getName()).log(Level.SEVERE, null, ex);
//            }
            System.out.println("raytracing.RayTracingCube.main()");
        }
    }

    private static void show() {

        form.getDrawable().add(new Lienzo.IDrawableInformer() {
            @Override
            public void paint(Graphics g) {
                g.drawImage(bufferedImage, 0, 0, form.getWidth(), form.getHeight(), form);
            }
        });
        form.setSize(bufferedImage.getHeight(), bufferedImage.getWidth());

        java.awt.EventQueue.invokeLater(new Runnable() {
            public void run() {

                form.setVisible(true);

            }
        });
    }

    private static BufferedImage generarImagen() {
        double time = System.currentTimeMillis();
//        double[] f = {1,2,3};
//        Sphere.VectorialProduct(new double[]{2,0,1}, new double[]{1,-1,3}, f);

        //double hhhh = Sphere.GetSphereIntersec(200, 200, 200, 100, 0, 0, 0, 1, 1, 1);
        double fMax = 400.0;
        BufferedImage bufferedImage = new BufferedImage(W, H, BufferedImage.TYPE_INT_RGB);
        Graphics g = bufferedImage.getGraphics();

        Color clrBackground = Color.black;
        g.setColor(clrBackground);
        g.fillRect(0, 0, W, H);
        Rectangle rect = new Rectangle(0, 0, W, H);

        ArrayList<Triangle3D> obj3dArrayList = new ArrayList<>();
        obj3dArrayList.add(new Triangle3D(100, 0, 0, 0, 100, 0, 0, 0, 0, 0, 0, 255));
        obj3dArrayList.add(new Triangle3D( 0, 100, 0, 100, 0, 0, 100, 100, 0,  0, 0, 255));
        Graphics graphics = g;

        // light position
        double lpx = (double) 200,
                lpy = (double) 200,
                lpz = (double) 200;
//        // light direction
//        double lvx = (double) -1,
//                lvy = (double) -1,
//                lvz = (double) -1;

        double d = 300;
//        double[] viewdirection = new double[]{1, 1, 0};

        //double[] u = {1,0,0};
        double[] Vup = {0, 0, 1};
        double[] COI = {0.0, 0.0, 0.0};
        double[] N = {px - COI[0], py - COI[1], pz - COI[2]};
        double N_mod = Triangle3D.modv(N);
        double[] n = {N[0] / N_mod, N[1] / N_mod, N[2] / N_mod};
        double[] U = {0, 0, 0};
        Triangle3D.VectorialProduct(Vup, n, U);
        double U_mod = Triangle3D.modv(U);
        double[] u = {U[0] / U_mod, U[1] / U_mod, U[2] / U_mod};
        double[] v = {0, 0, 0};
        Triangle3D.VectorialProduct(n, u, v);

        double[] vect = {0, 0, 0};

        for (int i = (int) rect.getX(); i < rect.getX() + rect.getWidth(); i++) {
            //double x = Sphere.GetCoord(rect.getX(), rect.getX() + rect.getWidth(), -fMax, fMax, i);
            for (int j = (int) rect.getY(); j < rect.getY() + rect.getHeight(); j++) {
                //double y = Sphere.GetCoord(rect.getY(), rect.getY() + rect.getHeight(), fMax, -fMax, j);
                double t = 1.0E10;
//                double vx = x - px, vy = y - py, vz = -pz;
                double vx = 0, vy = 0, vz = 0;

                //   Sphere.GetVector((int)((rect.getX() + rect.getWidth())/2), (int)((rect.getY() + rect.getHeight())/2), px, py, pz, n, rect, d, u, v, H, W, vect);
                Triangle3D.GetVector(i, j, px, py, pz, n, rect, d, u, v, H, W, vect);

                vx = vect[0];
                vy = vect[1];
                vz = vect[2];

                double mod_v = Triangle3D.modv(vx, vy, vz);
                vx = vx / mod_v;
                vy = vy / mod_v;
                vz = vz / mod_v;
                boolean bShadow = false;
                Triangle3D trianglehit = null;
                for (int k = 0; k < (int) obj3dArrayList.size(); k++) {
                    Triangle3D sphn = (Triangle3D) obj3dArrayList.get(k);
                    double taux = sphn.intersect(px, py, pz, vx, vy, vz);

//                            Sphere.GetSphereIntersec(sphn.cx, sphn.cy,
//                            sphn.cz,
//                            sphn.radius, px, py, pz, vx, vy, vz);
//                    if(taux!=-1)
//                        System.out.println(i+" "+j+" "+taux);
                    if (Double.isNaN(taux) || taux < 0) {
                        continue;
                    }
                    if (taux > 0 && taux < t) {
                        t = taux;
                        trianglehit = sphn;
                    }
                }

                Color color = new Color(10, 20, 10);
                if (trianglehit != null) {
                    double itx = px + t * vx, ity = py + t * vy, itz = pz
                            + t * vz;
                    // shadow
//                    double tauxla = Sphere.GetSphereIntersec(trianglehit.cx,
//                            trianglehit.cy, trianglehit.cz, trianglehit.radius,
//                            lpx, lpy, lpz, itx - lpx,
//                            ity - lpy, itz - lpz);

//                    boolean DARSKIDE = Sphere.GetSphereIntersecFar(trianglehit.cx,
//                            trianglehit.cy, trianglehit.cz, trianglehit.radius, itx,
//                            ity, itz, itx - lpx, ity - lpy, itz
//                            - lpz) > 0.0001;
//                    if (DARSKIDE == false) {
//                        for (int k = 0; k < (int) obj3dArrayList.size(); k++) {
//                            Triangle3D sphnb = (Sphere) (obj3dArrayList.get(k));
//                            if (sphnb != trianglehit) {
//                                double tauxlb = Sphere.GetSphereIntersec(sphnb.cx,
//                                        sphnb.cy, sphnb.cz, sphnb.radius, lpx,
//                                        lpy, lpz, itx - lpx, ity - lpy, itz
//                                        - lpz);
//                                if (tauxlb > 0 && tauxla < tauxlb) {
//                                    bShadow = true;
//                                    break;
//                                }
//                            } else {
//
//                            }
//                        }
//                    }
//                    double cost = Sphere.GetCosAngleV1V2(lvx, lvy, lvz, itx
//                            - trianglehit.cx, ity - trianglehit.cy, itz
//                            - trianglehit.cz);   //   direccion desde un punto
                    double cost = Triangle3D.GetCosAngleV1V2(lpx - itx, lpy - ity, lpz - itz, trianglehit.Nx, trianglehit.Ny, trianglehit.Nz);

                    if (cost < 0) {
                        cost = 0;
                    }
                    double fact = 1.0;
                    if (bShadow == true) {
                        fact = 0.5;
                    } else {
                        fact = 1.0;
                    }
                    double rgbR = trianglehit.clR * cost * fact;
                    double rgbG = trianglehit.clG * cost * fact;
                    double rgbB = trianglehit.clB * cost * fact;
                    color = new Color((int) rgbR, (int) rgbG, (int) rgbB);
                    //pen = new Pen(color);
                }
                //Brush brs = new SolidBrush(color);
                //graphics.FillRectangle(brs, i, j, 1, 1);
                //brs.Dispose();
                bufferedImage.setRGB(i, j, color.getRGB());

            }// for pixels lines
        }// for pixels columns

        System.out.println(System.currentTimeMillis() - time);
        return bufferedImage;

    }

    public static class Lienzo extends JFrame {

        @Override
        public void paint(Graphics g) {
            for (IDrawableInformer inst : drawable) {
                inst.paint(g);
            }
        }

        private ArrayList<IDrawableInformer> drawable = new ArrayList<>();

        private ArrayList<IDrawableInformer> getDrawable() {
            return drawable;
        }

        public interface IDrawableInformer {

            public void paint(Graphics g);
        }
    }

    public static class Triangle3D {

        public double p1x, p1y, p1z, p2x, p2y, p2z, p3x, p3y, p3z;
        private double Nx, Ny, Nz, d;
        private double edge0x, edge0y, edge0z, edge1x, edge1y, edge1z, edge2x, edge2y, edge2z;
        public double clR, clG, clB;

        public Triangle3D(double p1x, double p1y, double p1z, double p2x, double p2y, double p2z, double p3x, double p3y, double p3z, double clR, double clG, double clB) {
            this.p1x = p1x;
            this.p1y = p1y;
            this.p1z = p1z;
            this.p2x = p2x;
            this.p2y = p2y;
            this.p2z = p2z;
            this.p3x = p3x;
            this.p3y = p3y;
            this.p3z = p3z;
            this.clR = clR;
            this.clG = clG;
            this.clB = clB;

            //https://www.scratchapixel.com/lessons/3d-basic-rendering/ray-tracing-rendering-a-triangle/ray-triangle-intersection-geometric-solution
            //   falta calcular  el plano
            double[] ntemp = new double[3];
            Triangle3D.VectorialProduct(new double[]{p1x - p3x, p1y - p3y, p1z - p3z}, new double[]{p2x - p3x, p2y - p3y, p2z - p3z}, ntemp);
            this.Nx = ntemp[0];
            this.Ny = ntemp[1];
            this.Nz = ntemp[2];
            d = -(Nx * p1x + Ny * p1y + Nz * p1z);

            edge0x = p2x - p1x;
            edge0y = p2y - p1y;
            edge0z = p2z - p1z;

            edge1x = p3x - p2x;
            edge1y = p3y - p2y;
            edge1z = p3z - p2z;

            edge2x = p1x - p3x;
            edge2y = p1y - p3y;
            edge2z = p1z - p3z;
        }

        public static double t;
        public static double Px, Py, Pz;
        public static double[] vtemp1 = new double[3];
        public static double[] vtemp2 = new double[3];
        public static double[] vtemp3 = new double[3];
        public static double[] vtemp4 = new double[3];

        public double intersect(double px, double py, double pz, double vx, double vy, double vz) {

            t = -(Nx * px + Ny * py + Nz * pz + d) / (Nx * vx + Ny * vy + Nz * vz);
            Px = px + t * vx;
            Py = py + t * vy;
            Pz = pz + t * vz;

            vtemp4[0] = Nx;
            vtemp4[1] = Ny;
            vtemp4[2] = Nz;

            vtemp1[0] = edge0x;
            vtemp1[1] = edge0y;
            vtemp1[2] = edge0z;
            vtemp2[0] = Px - p1x;
            vtemp2[1] = Py - p1y;
            vtemp2[2] = Pz - p1z;
            Triangle3D.VectorialProduct(vtemp1, vtemp2, vtemp3);
            int s1 = (int) Math.signum(vtemp3[0] * vtemp4[0] + vtemp3[1] * vtemp4[1] + vtemp3[2] * vtemp4[2]);

            vtemp1[0] = edge1x;
            vtemp1[1] = edge1y;
            vtemp1[2] = edge1z;
            vtemp2[0] = Px - p2x;
            vtemp2[1] = Py - p2y;
            vtemp2[2] = Pz - p2z;
            Triangle3D.VectorialProduct(vtemp1, vtemp2, vtemp3);
            int s2 = (int) Math.signum(vtemp3[0] * vtemp4[0] + vtemp3[1] * vtemp4[1] + vtemp3[2] * vtemp4[2]);

            vtemp1[0] = edge2x;
            vtemp1[1] = edge2y;
            vtemp1[2] = edge2z;
            vtemp2[0] = Px - p3x;
            vtemp2[1] = Py - p3y;
            vtemp2[2] = Pz - p3z;
            Triangle3D.VectorialProduct(vtemp1, vtemp2, vtemp3);
            int s3 = (int) Math.signum(vtemp3[0] * vtemp4[0] + vtemp3[1] * vtemp4[1] + vtemp3[2] * vtemp4[2]);

            return (s1 == s2 && s2 == s3) ? t : Double.NaN;
        }

        public static double GetCosAngleV1V2(double v1x, double v1y, double v1z,
                double v2x, double v2y, double v2z) {

            return (v1x * v2x + v1y * v2y + v1z * v2z) / (modv(v1x, v1y, v1z)
                    * modv(v2x, v2y, v2z));
        }

        public static void VectorialProduct(double[] a, double[] b, double[] res) {

            res[0] = a[1] * b[2] - a[2] * b[1];
            res[1] = a[2] * b[0] - a[0] * b[2];
            res[2] = a[0] * b[1] - a[1] * b[0];
        }

        public static void GetVector(int i, int j, double px, double py, double pz, double[] n, Rectangle rect, double d, double[] u, double[] v, double H, double W, double[] vect) {

            // d=150;
            double X = rect.width;
            double Y = rect.height;
            double[] C = {px - n[0] * d, py - n[1] * d, pz - n[2] * d};
            double[] L = {C[0] - u[0] * W / 2 - v[0] * H / 2, C[1] - u[1] * W / 2 - v[1] * H / 2, C[2] - u[2] * W / 2 - v[2] * H / 2};

            vect[0] = L[0] + u[0] * i * W / X + v[0] * j * H / Y - px;
            vect[1] = L[1] + u[1] * i * W / X + v[1] * j * H / Y - py;
            vect[2] = L[2] + u[2] * i * W / X + v[2] * j * H / Y - pz;
        }

        public static double modv(double vx, double vy, double vz) {
            return Math.sqrt(vx * vx + vy * vy + vz * vz);
        }

        public static double modv(double[] v) {
            return Math.sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2]);
        }

    }

//    public static class Sphere {
//
////        private static double[] GetVector(int i, int j, double fMax, double[] d, double px, double py, double pz) {
////           double xx = ((p - i1) / (2*fMax)) * (2*fMax) - fMax;
////        }
////        public static double GetCoord(double i1, double i2, double w1, double w2,
////                double p) {
////            return ((p - i1) / (i2 - i1)) * (w2 - w1) + w1;
////        }
////        private static double[] GetVector(int i, int j, double fMax, double[] d, double px, double py, double pz, Rectangle rect) {
////            throw new UnsupportedOperationException("Not supported yet."); //To change body of generated methods, choose Tools | Templates.
////        }
////        private static double[] GetVector(int i, int j, double fMax, double[] d, double px, double py, double pz, Rectangle rect, double d0, double[] viewdirection) {
////            throw new UnsupportedOperationException("Not supported yet."); //To change body of generated methods, choose Tools | Templates.
////        }
//        public static void VectorialProduct(double[] a, double[] b, double[] res) {
//
//            res[0] = a[1] * b[2] - a[2] * b[1];
//            res[1] = a[2] * b[0] - a[0] * b[2];
//            res[2] = a[0] * b[1] - a[1] * b[0];
//        }
//
//        public static void GetVector(int i, int j, double px, double py, double pz, double[] n, Rectangle rect, double d, double[] u, double[] v, double H, double W, double[] vect) {
//
//            // d=150;
//            double X = rect.width;
//            double Y = rect.height;
//            double[] C = {px - n[0] * d, py - n[1] * d, pz - n[2] * d};
//            double[] L = {C[0] - u[0] * W / 2 - v[0] * H / 2, C[1] - u[1] * W / 2 - v[1] * H / 2, C[2] - u[2] * W / 2 - v[2] * H / 2};
//
//            vect[0] = L[0] + u[0] * i * W / X + v[0] * j * H / Y - px;
//            vect[1] = L[1] + u[1] * i * W / X + v[1] * j * H / Y - py;
//            vect[2] = L[2] + u[2] * i * W / X + v[2] * j * H / Y - pz;
//        }
//
//        public Sphere(double x, double y, double z, double r, double clr,
//                double clg, double clb) {
//            cx = x;
//            cy = y;
//            cz = z;
//            radius = r;
//            clR = clr;
//            clG = clg;
//            clB = clb;
//        }
//
//        public static double modv(double vx, double vy, double vz) {
//            return Math.sqrt(vx * vx + vy * vy + vz * vz);
//        }
//
//        public static double modv(double[] v) {
//            return Math.sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2]);
//        }
//
//        void Move(double vx, double vy, double vz) {
//            cx += vx;
//            cy += vy;
//            cz += vz;
//        }
//
//        void MoveTo(double vx, double vy, double vz) {
//            cx = vx;
//            cy = vy;
//            cz = vz;
//        }
//
//        void RotX(double angle) {
//            double y = cy * Math.cos(angle) - cz * Math.sin(angle);
//            double z = cy * Math.sin(angle) + cz * Math.cos(angle);
//            cy = y;
//            cz = z;
//        }
//
//        void RotY(double angle) {
//            double x = cx * Math.cos(angle) - cz * Math.sin(angle);
//            double z = cx * Math.sin(angle) + cz * Math.cos(angle);
//            cx = x;
//            cz = z;
//        }
//
//        public static double GetSphereIntersec(double cx, double cy, double cz,
//                double radius, double px, double py, double pz,
//                double vx, double vy, double vz) {
//            // x-xo 2 + y-yo 2 + z-zo 2 = r 2
//            // x,y,z = p+tv
//            // At2 + Bt + C = 0
//            double A = (vx * vx + vy * vy + vz * vz);
//            double B = 2.0 * (px * vx + py * vy + pz * vz - vx * cx - vy
//                    * cy - vz * cz);
//            double C = px * px - 2 * px * cx + cx * cx + py * py - 2 * py
//                    * cy + cy * cy + pz * pz - 2 * pz * cz + cz * cz
//                    - radius * radius;
//            double D = B * B - 4 * A * C;
//            double t = -1.0;
//            if (D >= 0) {
//                double t1 = (-B - Math.sqrt(D)) / (2.0 * A);
//                double t2 = (-B + Math.sqrt(D)) / (2.0 * A);
//                if (t1 < t2) {
//                    t = t1;
//                } else {
//                    t = t2;
//                }
//            }
//            return t;
//        }
//
//        public static double GetSphereIntersecFar(double cx, double cy, double cz,
//                double radius, double px, double py, double pz,
//                double vx, double vy, double vz) {
//            // x-xo 2 + y-yo 2 + z-zo 2 = r 2
//            // x,y,z = p+tv
//            // At2 + Bt + C = 0
//            double A = (vx * vx + vy * vy + vz * vz);
//            double B = 2.0 * (px * vx + py * vy + pz * vz - vx * cx - vy
//                    * cy - vz * cz);
//            double C = px * px - 2 * px * cx + cx * cx + py * py - 2 * py
//                    * cy + cy * cy + pz * pz - 2 * pz * cz + cz * cz
//                    - radius * radius;
//            double D = B * B - 4 * A * C;
//            double t = -1.0;
//            if (D >= 0) {
//                double t1 = (-B - Math.sqrt(D)) / (2.0 * A);
//                double t2 = (-B + Math.sqrt(D)) / (2.0 * A);
//                if (t1 > t2) {
//                    t = t1;
//                } else {
//                    t = t2;
//                }
//            }
//            return t;
//        }
//
//        public static double GetCosAngleV1V2(double v1x, double v1y, double v1z,
//                double v2x, double v2y, double v2z) {
//
//            return (v1x * v2x + v1y * v2y + v1z * v2z) / (modv(v1x, v1y, v1z)
//                    * modv(v2x, v2y, v2z));
//        }
//        public double cx, cy, cz, radius, clR, clG, clB;
//    }
}
