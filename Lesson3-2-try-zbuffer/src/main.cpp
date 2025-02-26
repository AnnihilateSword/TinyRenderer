#include <vector>
#include <cmath>
#include <cstdlib>
#include <limits>
#include "tgaimage.h"
#include "model.h"
#include "geometry.h"

const TGAColor white = TGAColor(255, 255, 255, 255);
const TGAColor red   = TGAColor(255, 0,   0,   255);
Model *model = NULL;
const int width  = 800;
const int height = 800;

void line(int x0, int y0, int x1, int y1, TGAImage &image, TGAColor color)
{
    bool steep = false;
    if (std::abs(x0-x1)<std::abs(y0-y1))
    {
        std::swap(x0, y0);
        std::swap(x1, y1);
        steep = true;
    }
    if (x0>x1)
    {
        std::swap(x0, x1);
        std::swap(y0, y1);
    }

    for (int x=x0; x<=x1; x++)
    {
        float t = (x-x0)/(float)(x1-x0);
        int y = y0*(1.-t) + y1*t;
        if (steep)
        {
            image.set(y, x, color);
        }
        else
        {
            image.set(x, y, color);
        }
    }
}

/**
 * @brief 计算质心坐标
 * @return 如果返回的坐标 (u,v,w) 满足 u,v,w >0，则点 P 在三角形内部
 */
Vec3f barycentric(Vec3f A, Vec3f B, Vec3f C, Vec3f P)
{
    Vec3f s[2];
    // 计算 [AC,AB,PA] 的 x 和 y 分量
    // s[0] 存储 x 分量的 AC,AB,PA
    // s[1] 存储 y 分量的 AC,AB,PA
    for (int i=2; i--; )
    {
        s[i][0] = C[i]-A[i]; // 向量 AC 的分量
        s[i][1] = B[i]-A[i]; // 向量 AB 的分量
        s[i][2] = A[i]-P[i]; // 向量 PA 的分量（注意负号）
    }

    // 这里 cross(s[0], s[1]) 相当于 AC 叉乘 AB，结果是 u 垂直于这两个向量
    // 它的 z 分量 u.z 表示三角形 ABC 所在平面的法向量的 z 分量
    // 叉积 |axb| = |a||b|sin(α)  值为 0 时角度为 0，即 ABC 三点共线
    // u.x = s[0].y · s[1].z - s[0].z · s[1].y
    // u.y = s[0].z · s[1].x - s[0].x · s[1].z
    // u.z = s[0].x · s[1].y - s[0].y · s[1].x
    Vec3f u = cross(s[0], s[1]);

    // 三点共线时，会导致 u[2] 为 0，此时返回 (-1,1,1)
    // 判断是否大于 0.01
    if (std::abs(u[2])>1e-2)
        //若 1-u-v，u，v全为大于 0 的数，表示点在三角形内部
        // P = (1-u-v)*A + u*B + v*C
        // u.x 和 u.y 是点 P 相对于三角形边 AB 和 AC 的贡献
        // 所以 u = u.y/u.z, v = u.x/u.z
        return Vec3f(1.f-(u.x+u.y)/u.z, u.y/u.z, u.x/u.z);

    // 在这种情况下生成负坐标，它将被栅格化器丢弃
    return Vec3f(-1,1,1);
}

void triangle(Vec3f *pts, float *zbuffer, TGAImage &image, TGAColor color)
{
    Vec2f bboxmin( std::numeric_limits<float>::max(),  std::numeric_limits<float>::max());
    Vec2f bboxmax(-std::numeric_limits<float>::max(), -std::numeric_limits<float>::max());
    // clamp 用于限制包围盒在图像范围内
    // clamp：表示图像的最大 x 和 y 坐标 (从 0 开始，所以需要减 1)。
    Vec2f clamp(image.get_width()-1, image.get_height()-1);
    for (int i=0; i<3; i++)
    {
        for (int j=0; j<2; j++)
        {
            // 计算三角形的包围盒
            // 保证不超出图像的左边界和上边界
            bboxmin[j] = std::max(0.f,      std::min(bboxmin[j], pts[i][j]));
            // 保证不超出图像的右边界和下边界
            bboxmax[j] = std::min(clamp[j], std::max(bboxmax[j], pts[i][j]));
        }
    }

    Vec3f P;

    // 遍历一个二维包围盒（Bounding Box）中的每一个像素点，并判断这些像素点是否位于给定的三角形内。
    // 如果像素点在三角形内，则计算其深度值（z 值），并进行深度测试（Z-Buffer 测试），
    // 最后更新帧缓冲区和深度缓冲区。

    // 遍历边框中的每一个像素点
    for (P.x=bboxmin.x; P.x<=bboxmax.x; P.x++)
    {
        for (P.y=bboxmin.y; P.y<=bboxmax.y; P.y++)
        {
            // 获取质心坐标
            Vec3f bc_screen  = barycentric(pts[0], pts[1], pts[2], P);
            // 质心坐标有一个负值，说明点在三角形外
            if (bc_screen.x<0 || bc_screen.y<0 || bc_screen.z<0) continue;
            P.z = 0;
            //计算 zbuffer，并且每个顶点的 z 值乘上对应的质心坐标分量
            // P.z = u.x*pts[0][2] + u.y*pts[1][2] + u.z*pts[2][2];
            // 其中 pts[0][2] 即 pts[0].z
            // 根据质心坐标 (u,v,w) 和三角形顶点的深度值 pts[i][2]，插值计算 P 的深度值：
            for (int i=0; i<3; i++)
                P.z += pts[i][2]*bc_screen[i];

            // 深度测试
            if (zbuffer[int(P.x+P.y*width)] < P.z)
            {
                zbuffer[int(P.x+P.y*width)] = P.z;
                image.set(P.x, P.y, color);
            }
        }
    }
}

Vec3f world2screen(Vec3f v)
{
    return Vec3f(int((v.x+1.)*width/2.+.5), int((v.y+1.)*height/2.+.5), v.z);
}

int main(int argc, char** argv)
{
    if (2==argc)
    {
        model = new Model(argv[1]);
    }
    else
    {
        model = new Model("../../res/obj/african_head/african_head.obj");
    }

    float *zbuffer = new float[width*height];
    for (int i=width*height; i--; zbuffer[i] = -std::numeric_limits<float>::max());

    TGAImage image(width, height, TGAImage::RGB);

    Vec3f light_dir(0,0,-1); // define light_dir

    for (int i=0; i<model->nfaces(); i++)
    {
        std::vector<int> face = model->face(i);
        Vec3f world_coords[3];
        for (int j=0; j<3; j++)
        {
            Vec3f v = model->vert(face[j]);
            world_coords[j]  = v;
        }
        Vec3f A = world_coords[2] - world_coords[0];
        Vec3f B = world_coords[1] - world_coords[0];
        // 叉积求面的法向量
        Vec3f n = Vec3f(A.y*B.z - A.z*B.y, A.z*B.x - A.x*B.z, A.x*B.y - A.y*B.x);

        n.normalize();
        // 计算光线强度
        float intensity = n*light_dir;
        if (intensity>0)
        {
            Vec3f pts[3];
            for (int i=0; i<3; i++) pts[i] = world2screen(model->vert(face[i]));
            triangle(pts, zbuffer, image, TGAColor(intensity*255, intensity*255, intensity*255, 255));
        }
    }

    image.flip_vertically(); // i want to have the origin at the left bottom corner of the image
    image.write_tga_file("output.tga");
    delete model;
    return 0;
}
