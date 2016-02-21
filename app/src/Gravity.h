#pragma once

#include <SFML/Graphics.hpp>

#include <vector>

class Gravity {
 public:
    Gravity();
    ~Gravity();

    void Add(int count);
    void Calc(double dt);
    void Draw(sf::RenderWindow &window);
    void GetCenter(float &x, float &y);

    int GetCount();

 private:
    struct Particle {
        double x, y;
        double vx, vy;
        double mass;
        double radius;
        double ax, ay;
    };

    double center_x;
    double center_y;
    double mass_all;

    int count;
    std::vector<Particle> particles;

};
