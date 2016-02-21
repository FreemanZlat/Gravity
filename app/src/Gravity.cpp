#include "Gravity.h"

#include <cmath>
#include <cstdlib>
#include <ctime>

static const double START_X = 640.0;
static const double START_Y = 360.0;
static const double START_SIZE = 400.0;

static const double GRAVITY_G = 6.673848e-11;
static const double GRAVITY_E = 0.1;
static const double GRAVITY_MASS = 100000000.0;
static const double GRAVITY_MASS_BIG = 2000000000000.0;
static const double MAX_DISTANCE = 4096;

static double Random(double q)
{
    return q * static_cast<double>(rand()) / (static_cast<double>(RAND_MAX) + 1.0);
}

Gravity::Gravity() :
        center_x(0.0),
        center_y(0.0),
        mass_all(0.0),
        count(1)
{
    srand((unsigned) time(nullptr));

    Particle p;
    p.x = START_X;
    p.y = START_Y;
    p.vx = 0.0;
    p.vy = 0.0;
    p.mass = GRAVITY_MASS_BIG;
    p.radius = 0.007 * sqrt(sqrt(p.mass));
    this->particles.push_back(p);
}

Gravity::~Gravity()
{
}

void Gravity::Add(int count)
{
    for (int i = 0; i < count; ++i)
    {
        Particle p;

        double x = Random(2.0 * START_SIZE) - START_SIZE;
        double y = Random(2.0 * START_SIZE) - START_SIZE;
        double r = sqrt(x * x + y * y);

        p.x = START_X + x;
        p.y = START_Y + y;

        p.vx = y * 0.02 / sqrt(r);
        p.vy = -x * 0.02 / sqrt(r);

        p.mass = GRAVITY_MASS;
        p.radius = 0.007 * sqrt(sqrt(p.mass));

        particles.push_back(p);
    }

    this->count = this->particles.size();
}

void Gravity::Calc(double dt)
{
    for (auto &p : this->particles)
    {
        p.ax = 0.0;
        p.ay = 0.0;
    }

    this->center_x = 0.0;
    this->center_y = 0.0;
    this->mass_all = 0.0;

    for (int i = 0; i < this->count; ++i)
    {
        Particle &p1 = this->particles[i];

        for (int j = i + 1; j < this->count; ++j)
        {
            Particle &p2 = this->particles[j];

            double dx = p2.x - p1.x;
            double dy = p2.y - p1.y;
            double r2 = (dx * dx) + (dy * dy);
            double min_r = p1.radius + p2.radius;

            if (r2 <= min_r * min_r)
            {
                p1.x = p1.x * p1.mass + p2.x * p2.mass;
                p1.y = p1.y * p1.mass + p2.y * p2.mass;
                p1.vx = p1.vx * p1.mass + p2.vx * p2.mass;
                p1.vy = p1.vy * p1.mass + p2.vy * p2.mass;

                p1.mass += p2.mass;
                p1.radius = 0.007 * sqrt(sqrt(p1.mass));

                p1.x /= p1.mass;
                p1.y /= p1.mass;
                p1.vx /= p1.mass;
                p1.vy /= p1.mass;

                this->count--;
                p2 = this->particles[this->count];
                this->particles.erase(this->particles.begin() + this->count);

                j--;
                continue;
            }

            double k = 1.0 / (r2 * sqrt(r2) + GRAVITY_E);

            p1.ax += p2.mass * k * dx;
            p1.ay += p2.mass * k * dy;

            p2.ax -= p1.mass * k * dx;
            p2.ay -= p1.mass * k * dy;
        }

        p1.vx += GRAVITY_G * p1.ax * dt;
        p1.vy += GRAVITY_G * p1.ay * dt;

        p1.x += p1.vx * dt;
        p1.y += p1.vy * dt;

        this->center_x += p1.x * p1.mass;
        this->center_y += p1.y * p1.mass;
        this->mass_all += p1.mass;
    }

    this->center_x /= this->mass_all;
    this->center_y /= this->mass_all;

    for (int i = 0; i < this->count; ++i)
    {
        double dx = this->particles[i].x - this->center_x;
        double dy = this->particles[i].y - this->center_y;

        if (dx * dx + dy * dy > MAX_DISTANCE * MAX_DISTANCE)
        {
            this->count--;
            this->particles[i] = this->particles[this->count];
            this->particles.erase(this->particles.begin() + this->count);
            i--;
        }
    }
}

void Gravity::Draw(sf::RenderWindow &window)
{
    sf::CircleShape circle(1, 16);
    for (auto &p : this->particles)
    {
        circle.setPosition(p.x - p.radius, p.y - p.radius);
        circle.setRadius(p.radius);
        window.draw(circle);
    }
}

void Gravity::GetCenter(float &x, float &y)
{
    x = this->center_x;
    y = this->center_y;
}

int Gravity::GetCount()
{
    return this->count;
}
