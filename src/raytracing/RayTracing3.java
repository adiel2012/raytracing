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
public class RayTracing3 {

    public static int W = 400, H = 400;
    public static Lienzo form = new Lienzo();
    public static double[] COI = {0.0, 0.0, 90.0};
    // viewer position
    public static double px = (double) 0,
            py = (double) 0,
            pz = (double) 500;
    /**
     * @param args the command line arguments
     */
    //   http://www.codeproject.com/Articles/19732/Simple-Ray-Tracing-in-C

    public static BufferedImage bufferedImage = null;

    public static void main(String[] args) {
        bufferedImage = generarImagen();
        show();
        double angle = Math.PI / 360;
        while (true) {
            bufferedImage = generarImagen();
            form.repaint();
            // pz -= 10 ;

            double x = px - COI[0], y = py - COI[1], z = pz - COI[2];

            double xm = x * Math.cos(angle) - z * Math.sin(angle);
            double zm = z * Math.cos(angle) + x * Math.sin(angle);

            px = xm + COI[0];
            py = y + COI[1];
            pz = zm + COI[2];

            System.out.println("raytracing.RayTracing2.main()");
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

        ArrayList<Sphere> obj3dArrayList = new ArrayList<>();
        obj3dArrayList.add(new RayTracing3.Sphere(0.0, 0.0, 90.0, 100.0, 0.0, 0.0, 255.0));
        obj3dArrayList.add(new RayTracing3.Sphere(-160.0, -120.0, -110.0, 15.0, 255.0, 0.0, 0.0));
        obj3dArrayList.add(new RayTracing3.Sphere(-140.0, -140.0, -150.0, 20.0, 255.0, 200.0, 0.0));
        Graphics graphics = g;

        // light position
        double lpx = (double) 200,
                lpy = (double) 200,
                lpz = (double) 200;
        // light direction
        double lvx = (double) -1,
                lvy = (double) -1,
                lvz = (double) -1;

        double d = 100;
//        double[] viewdirection = new double[]{1, 1, 0};

        //double[] u = {1,0,0};
        double[] Vup = {-1, 1, 0};

        double[] N = {px - COI[0], py - COI[1], pz - COI[2]};
        double N_mod = Sphere.modv(N);
        double[] n = {N[0] / N_mod, N[1] / N_mod, N[2] / N_mod};
        double[] U = {0, 0, 0};
        Sphere.VectorialProduct(Vup, n, U);
        double U_mod = Sphere.modv(U);
        double[] u = {U[0] / U_mod, U[1] / U_mod, U[2] / U_mod};
        double[] v = {0, 0, 0};
        Sphere.VectorialProduct(n, u, v);

        double[] vect = {0, 0, 0};

        for (int i = (int) rect.getX(); i < rect.getX() + rect.getWidth(); i++) {
            double x = Sphere.GetCoord(rect.getX(), rect.getX() + rect.getWidth(), -fMax, fMax, i);
            for (int j = (int) rect.getY(); j < rect.getY() + rect.getHeight(); j++) {
                double y = Sphere.GetCoord(rect.getY(), rect.getY() + rect.getHeight(), fMax, -fMax, j);
                double t = 1.0E10;
                double vx = x - px, vy = y - py, vz = -pz;

                //   Sphere.GetVector((int)((rect.getX() + rect.getWidth())/2), (int)((rect.getY() + rect.getHeight())/2), px, py, pz, n, rect, d, u, v, H, W, vect);
                Sphere.GetVector(i, j, px, py, pz, n, rect, d, u, v, H, W, vect);

                vx = vect[0];
                vy = vect[1];
                vz = vect[2];

                double mod_v = Sphere.modv(vx, vy, vz);
                vx = vx / mod_v;
                vy = vy / mod_v;
                vz = vz / mod_v;
                boolean bShadow = false;
                Sphere spherehit = null;
                for (int k = 0; k < (int) obj3dArrayList.size(); k++) {
                    Sphere sphn = (Sphere) obj3dArrayList.get(k);
                    double taux = Sphere.GetSphereIntersec(sphn.cx, sphn.cy,
                            sphn.cz,
                            sphn.radius, px, py, pz, vx, vy, vz);
//                    if(taux!=-1)
//                        System.out.println(i+" "+j+" "+taux);
                    if (taux < 0) {
                        continue;
                    }
                    if (taux > 0 && taux < t) {
                        t = taux;
                        spherehit = sphn;
                    }
                }
                Color color = new Color(10, 20, 10);
                if (spherehit != null) {
                    double itx = px + t * vx, ity = py + t * vy, itz = pz
                            + t * vz;
                    // shadow
                    double tauxla = Sphere.GetSphereIntersec(spherehit.cx,
                            spherehit.cy, spherehit.cz, spherehit.radius,
                            lpx, lpy, lpz, itx - lpx,
                            ity - lpy, itz - lpz);
                    boolean DARSKIDE = RayTracing2.Sphere.GetSphereIntersecFar(spherehit.cx,
                            spherehit.cy, spherehit.cz, spherehit.radius, itx,
                            ity, itz, itx - lpx, ity - lpy, itz
                            - lpz) > 0.0001;
                    if (DARSKIDE == false) {
                        for (int k = 0; k < (int) obj3dArrayList.size(); k++) {
                            Sphere sphnb = (Sphere) (obj3dArrayList.get(k));
                            if (sphnb != spherehit) {
                                double tauxlb = Sphere.GetSphereIntersec(sphnb.cx,
                                        sphnb.cy, sphnb.cz, sphnb.radius, lpx,
                                        lpy, lpz, itx - lpx, ity - lpy, itz
                                        - lpz);
                                if (tauxlb > 0 && tauxla < tauxlb) {
                                    bShadow = true;
                                    break;
                                }
                            }
                        }
                    }
                    double cost = Sphere.GetCosAngleV1V2(lvx, lvy, lvz, itx
                            - spherehit.cx, ity - spherehit.cy, itz
                            - spherehit.cz);
                    if (cost < 0) {
                        cost = 0;
                    }
                    double fact = 1.0;
                    if (bShadow == true) {
                        fact = 0.5;
                    } else {
                        fact = 1.0;
                    }
                    double rgbR = spherehit.clR * cost * fact;
                    double rgbG = spherehit.clG * cost * fact;
                    double rgbB = spherehit.clB * cost * fact;
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

    public static class Sphere {

//        private static double[] GetVector(int i, int j, double fMax, double[] d, double px, double py, double pz) {
//           double xx = ((p - i1) / (2*fMax)) * (2*fMax) - fMax;
//        }
        public static double GetCoord(double i1, double i2, double w1, double w2,
                double p) {
            return ((p - i1) / (i2 - i1)) * (w2 - w1) + w1;
        }

//        private static double[] GetVector(int i, int j, double fMax, double[] d, double px, double py, double pz, Rectangle rect) {
//            throw new UnsupportedOperationException("Not supported yet."); //To change body of generated methods, choose Tools | Templates.
//        }
//        private static double[] GetVector(int i, int j, double fMax, double[] d, double px, double py, double pz, Rectangle rect, double d0, double[] viewdirection) {
//            throw new UnsupportedOperationException("Not supported yet."); //To change body of generated methods, choose Tools | Templates.
//        }
        private static void VectorialProduct(double[] a, double[] b, double[] res) {

            res[0] = a[1] * b[2] - a[2] * b[1];
            res[1] = a[2] * b[0] - a[0] * b[2];
            res[2] = a[0] * b[1] - a[1] * b[0];
        }

        private static void GetVector(int i, int j, double px, double py, double pz, double[] n, Rectangle rect, double d, double[] u, double[] v, double H, double W, double[] vect) {

            // d=150;
            double X = rect.width;
            double Y = rect.height;
            double[] C = {px - n[0] * d, py - n[1] * d, pz - n[2] * d};
            double[] L = {C[0] - u[0] * W / 2 - v[0] * H / 2, C[1] - u[1] * W / 2 - v[1] * H / 2, C[2] - u[2] * W / 2 - v[2] * H / 2};

            vect[0] = L[0] + u[0] * i * W / X + v[0] * j * H / Y - px;
            vect[1] = L[1] + u[1] * i * W / X + v[1] * j * H / Y - py;
            vect[2] = L[2] + u[2] * i * W / X + v[2] * j * H / Y - pz;
        }

        public Sphere(double x, double y, double z, double r, double clr,
                double clg, double clb) {
            cx = x;
            cy = y;
            cz = z;
            radius = r;
            clR = clr;
            clG = clg;
            clB = clb;
        }

        public static double modv(double vx, double vy, double vz) {
            return Math.sqrt(vx * vx + vy * vy + vz * vz);
        }

        public static double modv(double[] v) {
            return Math.sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2]);
        }

        void Move(double vx, double vy, double vz) {
            cx += vx;
            cy += vy;
            cz += vz;
        }

        void MoveTo(double vx, double vy, double vz) {
            cx = vx;
            cy = vy;
            cz = vz;
        }

        void RotX(double angle) {
            double y = cy * Math.cos(angle) - cz * Math.sin(angle);
            double z = cy * Math.sin(angle) + cz * Math.cos(angle);
            cy = y;
            cz = z;
        }

        void RotY(double angle) {
            double x = cx * Math.cos(angle) - cz * Math.sin(angle);
            double z = cx * Math.sin(angle) + cz * Math.cos(angle);
            cx = x;
            cz = z;
        }

        public static double GetSphereIntersec(double cx, double cy, double cz,
                double radius, double px, double py, double pz,
                double vx, double vy, double vz) {
            // x-xo 2 + y-yo 2 + z-zo 2 = r 2
            // x,y,z = p+tv
            // At2 + Bt + C = 0
            double A = (vx * vx + vy * vy + vz * vz);
            double B = 2.0 * (px * vx + py * vy + pz * vz - vx * cx - vy
                    * cy - vz * cz);
            double C = px * px - 2 * px * cx + cx * cx + py * py - 2 * py
                    * cy + cy * cy + pz * pz - 2 * pz * cz + cz * cz
                    - radius * radius;
            double D = B * B - 4 * A * C;
            double t = -1.0;
            if (D >= 0) {
                double t1 = (-B - Math.sqrt(D)) / (2.0 * A);
                double t2 = (-B + Math.sqrt(D)) / (2.0 * A);
                if (t1 < t2) {
                    t = t1;
                } else {
                    t = t2;
                }
            }
            return t;
        }

        public static double GetSphereIntersecFar(double cx, double cy, double cz,
                double radius, double px, double py, double pz,
                double vx, double vy, double vz) {
            // x-xo 2 + y-yo 2 + z-zo 2 = r 2
            // x,y,z = p+tv
            // At2 + Bt + C = 0
            double A = (vx * vx + vy * vy + vz * vz);
            double B = 2.0 * (px * vx + py * vy + pz * vz - vx * cx - vy
                    * cy - vz * cz);
            double C = px * px - 2 * px * cx + cx * cx + py * py - 2 * py
                    * cy + cy * cy + pz * pz - 2 * pz * cz + cz * cz
                    - radius * radius;
            double D = B * B - 4 * A * C;
            double t = -1.0;
            if (D >= 0) {
                double t1 = (-B - Math.sqrt(D)) / (2.0 * A);
                double t2 = (-B + Math.sqrt(D)) / (2.0 * A);
                if (t1 > t2) {
                    t = t1;
                } else {
                    t = t2;
                }
            }
            return t;
        }

        public static double GetCosAngleV1V2(double v1x, double v1y, double v1z,
                double v2x, double v2y, double v2z) {
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
            return (v1x * v2x + v1y * v2y + v1z * v2z) / (modv(v1x, v1y, v1z)
                    * modv(v2x, v2y, v2z));
        }
        public double cx, cy, cz, radius, clR, clG, clB;
    }

}
