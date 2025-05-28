#ifndef CAMERA_H
#define CAMERA_H



#include "hittable.h"
#include "material.h"

#pragma once
#include "rtweekend.h"
#include "hittable.h"
#include "color.h"
#include <vector>
#include <thread>
#include <atomic>
#include <mutex>
#include <iostream>

double halton(int index, int base)
{
    double result = 0.0;
    double f = 1.0;
    while (index > 0)
    {
        f /= base;
        result += f * (index % base);
        index /= base;
    }
    return result;
}
class camera
{
public:
    /* Public Camera Parameters Here */
    double aspect_ratio = 1.0;
    int image_width = 100;
    int samples_per_pixel = 10;
    int max_depth = 50;

    double vfov = 90;
    point3 lookfrom = point3(0, 0, 0); // where we at
    point3 lookat = point3(0, 0, -1);  // what we see ;
    vec3 vup = vec3(0, 1, 0);          //

    double defocus_angle = 0;
    double focus_dist = 10;

    void render(const hittable &world)
    {
        initialize();

        std::cout << "P3\n"
                  << image_width << " " << image_height << "\n255\n";

        for (int j = 0; j < image_height; j++)
        {
            std::clog << "\rScanlines remaining: " << (image_height - j) << ' ' << std::flush;
            for (int i = 0; i < image_width; i++)
            {
                color pixel_color(0, 0, 0);
                for (int sample = 0; sample < samples_per_pixel; sample++)
                {
                    ray r = get_ray(i, j);
                    pixel_color += ray_color(r,max_depth, world);
                }
                write_color(std::cout, pixel_samples_scale * pixel_color);
            }
        }

        std::clog << "\rDone.                 \n";
    }

private:
    /* Private Camera Variables Here */
    int image_height;
    double pixel_samples_scale;
    point3 center;
    point3 pixel00_loc;
    vec3 pixel_delta_u;
    vec3 pixel_delta_v;

    vec3 u, v, w; // Camera frame basis vectors

    vec3 defocus_disk_u;
    vec3 defocus_disk_v;

    void initialize()
    {
        image_height = int(image_width / aspect_ratio);
        image_height = (image_height < 1) ? 1 : image_height;

        pixel_samples_scale = 1.0 / samples_per_pixel;

        center = lookfrom;

        auto theta = degrees_to_radians(vfov);
        auto h = std::tan(theta / 2);

        auto viewport_height = 2 * h * focus_dist; // focal_length is -z .
        auto viewport_width = viewport_height * (double(image_width) / image_height);

        // Caculate the u,v,w unit basis vectors  for th e camera coordinate frame.

        w = unit_vector(lookfrom - lookat);
        u = unit_vector(cross(vup, w));
        v = cross(w, u);

        // Calculate the vectors across the horizontal and down the vertical viewport edges.
        auto viewport_u = viewport_width * u;
        auto viewport_v = viewport_height * -v;

        // Calculate the horizontal and vertical delta vectors from pixel to pixel.

        pixel_delta_u = viewport_u / image_width;
        pixel_delta_v = viewport_v / image_height;

        // Calculate the location of the upper left pixel.
        auto viewport_upper_left = center - (focus_dist * w) - viewport_u / 2 - viewport_v / 2;
        pixel00_loc = viewport_upper_left + 0.5 * (pixel_delta_u + pixel_delta_v);

        // Calculate the camera defocus disk basis vectors;
        auto defocus_radius = focus_dist * std::tan(degrees_to_radians(defocus_angle / 2));
        // defocus_angle 有点像 vfov,然后这个defocus_radius 类似之前vfov的h，就是那个高度。
        defocus_disk_u = u * defocus_radius;
        defocus_disk_v = v * defocus_radius;
    }

    ray get_ray(int i, int j) const
    {
        // Construct a camera ray originating from the defocus disk and directed at a randomly
        // sampled point around the pixel location i, j.
        auto offset = sample_square();
        auto pixel_sample = pixel00_loc + ((i + offset.x()) * pixel_delta_u) + ((j + offset.y()) * pixel_delta_v);

        // add defocus in ray_oringin ,相当于一个从0的disk变成有radius的一个defocus disk。
        auto ray_origin = (defocus_angle <= 0) ? center : defocus_disk_sample(); // 从圆盘中随机挑一个点
        auto ray_direction = pixel_sample - ray_origin;
        return ray(ray_origin, ray_direction);
    }

    vec3 sample_square() const
    {
        // 使用Halton序列（需全局维护采样计数器）
        static int sample_count = 0;
        double x = halton(sample_count, 2) - 0.5; // Base 2
        double y = halton(sample_count, 3) - 0.5; // Base 3
        sample_count++;
        return vec3(x, y, 0);
        // return vec3(random_double() - 0.5, random_double() - 0.5, 0);
    }
    point3 defocus_disk_sample() const
    {
        // 返回一个相机defocus圆盘
        auto p = random_in_unit_disk();
        return center + (p[0] * defocus_disk_u) + (p[1] * defocus_disk_v);
        //如果用center +p*radius，可能会因为相机倾斜而导致生成点不正确，所以这里用的是基向量，就能保证正交。
    }
    color ray_color(const ray &r, int depth, const hittable &world)
    {
        if (depth <= 0)
        {
            return color(0, 0, 0);
        }

        hit_record rec;

        if (world.hit(r, interval(0.001, infinity), rec))
        {
            ray scattered;
            color attenuation;
            if (rec.mat->scatter(r, rec, attenuation, scattered))
            {
                return attenuation * ray_color(scattered, depth - 1, world);
            }
            return color(0, 0, 0);
        }

        vec3 unit_direction = unit_vector(r.direction());
        auto a = 0.5 * (unit_direction.y() + 1.0);
        return (1.0 - a) * color(1.0, 1.0, 1.0) + a * color(0.5, 0.7, 1.0);
    }
};

#endif