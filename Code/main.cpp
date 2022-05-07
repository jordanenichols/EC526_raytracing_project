#include "rtweekend.h"

#include "camera.h"
#include "color.h"
#include "hittable_list.h"
#include "material.h"
#include "sphere.h"

#include <iostream>
#include <fstream>
#include <time.h>
#include <omp.h>
#define PRINT_TIME 0

using namespace std;

color ray_color(const ray &r, const hittable &world, int depth)
{
    hit_record rec;

    // If we've exceeded the ray bounce limit, no more light is gathered.
    if (depth <= 0)
        return color(0, 0, 0);

    if (world.hit(r, 0.001, infinity, rec))
    {
        ray scattered;
        color attenuation;
        if (rec.mat_ptr->scatter(r, rec, attenuation, scattered))
            return attenuation * ray_color(scattered, world, depth - 1);
        return color(0, 0, 0);
    }

    vec3 unit_direction = unit_vector(r.direction());
    auto t = 0.5 * (unit_direction.y() + 1.0);
    return (1.0 - t) * color(1.0, 1.0, 1.0) + t * color(0.5, 0.7, 1.0);
}
hittable_list my_scene()
{
    hittable_list world;

    auto ground_material = make_shared<lambertian>(color(0.5, 0.5, 0.5));
    world.add(make_shared<sphere>(point3(0, -1000, 0), 1000, ground_material));

    auto material1 = make_shared<dielectric>(1.5);
    world.add(make_shared<sphere>(point3(0, 1, 0), 1.0, material1));

    auto material2 = make_shared<lambertian>(color(0.4, 0.2, 0.1));
    world.add(make_shared<sphere>(point3(-4, 1, 0), 1.0, material2));

    auto material3 = make_shared<metal>(color(0.7, 0.6, 0.5), 0.0);
    world.add(make_shared<sphere>(point3(4, 1, 0), 1.0, material3));

    return world;
}

hittable_list random_scene()
{
    hittable_list world;

    auto ground_material = make_shared<lambertian>(color(0.5, 0.5, 0.5));
    world.add(make_shared<sphere>(point3(0, -1000, 0), 1000, ground_material));

    for (int a = -11; a < 11; a++)
    {
        for (int b = -11; b < 11; b++)
        {
            auto choose_mat = random_double();
            point3 center(a + 0.9 * random_double(), 0.2, b + 0.9 * random_double());

            if ((center - point3(4, 0.2, 0)).length() > 0.9)
            {
                shared_ptr<material> sphere_material;

                if (choose_mat < 0.8)
                {
                    // diffuse
                    auto albedo = color::random() * color::random();
                    sphere_material = make_shared<lambertian>(albedo);
                    world.add(make_shared<sphere>(center, 0.2, sphere_material));
                }
                else if (choose_mat < 0.95)
                {
                    // metal
                    auto albedo = color::random(0.5, 1);
                    auto fuzz = random_double(0, 0.5);
                    sphere_material = make_shared<metal>(albedo, fuzz);
                    world.add(make_shared<sphere>(center, 0.2, sphere_material));
                }
                else
                {
                    // glass
                    sphere_material = make_shared<dielectric>(1.5);
                    world.add(make_shared<sphere>(center, 0.2, sphere_material));
                }
            }
        }
    }

    auto material1 = make_shared<dielectric>(1.5);
    world.add(make_shared<sphere>(point3(0, 1, 0), 1.0, material1));

    auto material2 = make_shared<lambertian>(color(0.4, 0.2, 0.1));
    world.add(make_shared<sphere>(point3(-4, 1, 0), 1.0, material2));

    auto material3 = make_shared<metal>(color(0.7, 0.6, 0.5), 0.0);
    world.add(make_shared<sphere>(point3(4, 1, 0), 1.0, material3));

    return world;
}

void write_buf_to_ppm(string filename, color **buf, const int image_width, const int image_height) {
    color c;
    ofstream outfile;
    outfile.open(filename);
    outfile << "P3\n"
                << image_width << ' ' << image_height << "\n255\n";
    for(int j = image_height - 1; j >= 0; --j) {
        for(int i = 0; i < image_width; i++) {
            c = buf[i][j];
            outfile << c[0] << ' ' << c[1] << ' ' << c[2] << '\n';
        }
    }
    outfile.close();
}

void generate_image_scalar(string filename, hittable_list world, camera cam, float aspect_ratio,
                           const int image_width, const int image_height, int samples_per_pixel, int max_depth)
{
    clock_t time;
    color** buf = new color*[image_width];
    for(int i = 0; i < image_width; ++i)
        buf[i] = new color[image_height];

    time = -clock();
    for (int j = image_height - 1; j >= 0; --j)
    {
        for (int i = 0; i < image_width; ++i)
        {
            color pixel_color(0, 0, 0);
            for (int s = 0; s < samples_per_pixel; ++s)
            {
                auto u = (i + random_double()) / (image_width - 1);
                auto v = (j + random_double()) / (image_height - 1);
                ray r = cam.get_ray(u, v);
                pixel_color += ray_color(r, world, max_depth);
            }
            write_color_to_buf(buf, i, j, pixel_color, samples_per_pixel);
        }
    }
    time += clock();
    write_buf_to_ppm(filename, buf, image_width, image_height);
    if (PRINT_TIME)
        std::cout << "Scalar:\t" << time << std::endl;
}

void generate_image_omp(string filename, hittable_list world, camera cam, float aspect_ratio,
                           const int image_width, const int image_height, int samples_per_pixel, int max_depth)
{
    clock_t time;
    color** buf = new color*[image_width];
    for(int i = 0; i < image_width; ++i)
        buf[i] = new color[image_height];

    time = -clock();
    #pragma omp parallel for
    for (int j = image_height - 1; j >= 0; --j)
    {
        for (int i = 0; i < image_width; ++i)
        {
            color pixel_color(0, 0, 0);
            for (int s = 0; s < samples_per_pixel; ++s)
            {
                auto u = (i + random_double()) / (image_width - 1);
                auto v = (j + random_double()) / (image_height - 1);
                ray r = cam.get_ray(u, v);
                pixel_color += ray_color(r, world, max_depth);
            }
            write_color_to_buf(buf, i, j, pixel_color, samples_per_pixel);
        }
    }
    time += clock();

    write_buf_to_ppm(filename, buf, image_width, image_height);
    if(PRINT_TIME)
        std::cout << "OpenMP:\t" << time << std::endl;
}

int main()
{
    // image values
    const auto aspect_ratio = 16.0 / 9.0;
    const int image_width = 640;
    const int image_height = static_cast<int>(image_width / aspect_ratio);
    const int samples_per_pixel = 10;
    const int max_depth = 10;

    // camera values
    point3 lookfrom(13, 2, 3);
    point3 lookat(0, 0, 0);
    vec3 vup(0, 1, 0);
    auto dist_to_focus = 10.0;
    auto aperture = 0.1;
    camera cam(lookfrom, lookat, vup, 20, aspect_ratio, aperture, dist_to_focus);


    auto world = random_scene();

    generate_image_scalar("test_scal.ppm", world, cam, aspect_ratio, image_width, image_height, samples_per_pixel, max_depth);
    generate_image_omp("test_omp.ppm", world, cam, aspect_ratio, image_width, image_height, samples_per_pixel, max_depth);
    return 0;
}
