#pragma once

#include <SFML/Graphics.hpp>
#include <glm/glm.hpp>

#include <vector>

class Gravity
{
 public:
    Gravity();
    ~Gravity();

    typedef glm::dvec2 vec;

    void Add(int count);
    void DoGravity(double dt);
    void DoCollisions();
    void DoClean();
    void Draw(sf::RenderWindow &window);
    void GetCenter(float &x, float &y);

    int GetCount();

 private:
    struct Info
    {
        double mass;
        double radius;
    };

    struct Particles
    {
        std::vector<vec> pos;
        std::vector<vec> vel;
        std::vector<Info> info;

        int count;
        vec center;
        double mass_all;

        void add(vec pos, vec vel, double mass);
        void del(int idx);
    };

    Particles particles;

    void calc_gravity(std::vector<vec> &pos, std::vector<vec> &accel);
    void calc_gravity_euler1(double dt);
    void calc_gravity_euler2(double dt);

    static double Random(double q);
};
