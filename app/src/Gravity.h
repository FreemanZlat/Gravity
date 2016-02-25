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
    void CalcGravity(double dt);
    void CalcCollisions();
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
        std::vector<vec> accel;
        std::vector<Info> info;

        int count;
        vec center;
        double mass_all;

        void add(vec pos, vec vel, double mass);
        void del(int idx);
    };

    Particles particles;

    static double Random(double q);
};
