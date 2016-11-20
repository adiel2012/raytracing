/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package reaytracing;

import java.awt.Color;
import java.awt.Graphics;
import java.awt.Rectangle;
import java.awt.image.BufferedImage;
import java.util.ArrayList;
import javax.swing.JFrame;
import trash.Lienzo;

/**
 *
 * @author adiel
 */
public class ReayTracing {

    /**
     * @param args the command line arguments
     */
    //   http://www.codeproject.com/Articles/19732/Simple-Ray-Tracing-in-C
    public static void main(String[] args) {
        int W = 200, H = 200;
        BufferedImage bufferedImage = new BufferedImage(W, H, BufferedImage.TYPE_INT_RGB);
        Graphics g = bufferedImage.getGraphics();

        Color clrBackground = Color.black;
        g.setColor(clrBackground);
        g.fillRect(0, 0, W, H);
        Rectangle rect = new Rectangle(0, 0, 200, 200);

        ArrayList<Sphere> obj3dArrayList = new ArrayList<>();
        obj3dArrayList.add(new Sphere(0.0, 0.0, 90.0, 100.0, 0.0, 0.0, 255.0));
        obj3dArrayList.add(new Sphere(-180.0, -130.0, -110.0, 15.0, 255.0, 0.0, 0.0));
        obj3dArrayList.add(new Sphere(-140.0, -140.0, -150.0, 20.0, 255.0, 200.0, 0.0));
        Graphics graphics = g;

        // viewer position
        double px = (double) 0,
                py = (double) 0,
                pz = (double) 500;
        // light position
        double lpx = (double) 200,
                lpy = (double) 200,
                lpz = (double) 200;
        // light direction
        double lvx = (double) -1,
                lvy = (double) -1,
                lvz = (double) -1;
        double fMax = 200.0;

        for (int i = (int) rect.getX(); i < rect.getX() + rect.getWidth(); i++) {
            double x = Sphere.GetCoord(rect.getX(), rect.getX()+rect.getWidth(), -fMax, fMax, i);
            for (int j = (int) rect.getY(); j < rect.getY()+rect.getHeight(); j++) {
                double y = Sphere.GetCoord(rect.getY(), rect.getY() + rect.getHeight(), fMax, -fMax, j);
                double t = 1.0E10;
                double vx = x - px, vy = y - py, vz = -pz;
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

        show(bufferedImage);
    }

    private static void show(BufferedImage bufferedImage) {
        Lienzo form = new Lienzo();
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

}
