#include "Gravity.h"

#include <SFML/Graphics.hpp>

#include <cstdio>
#include <sstream>

int main()
{
    // Fixing Eclipse x64 bug
    setvbuf(stdout, NULL, _IONBF, 0);
    setvbuf(stderr, NULL, _IONBF, 0);

    sf::Font font;
    if (!font.loadFromFile("../data/sansation.ttf"))
        return EXIT_FAILURE;

    sf::ContextSettings settings;
    settings.antialiasingLevel = 2;

    sf::RenderWindow window(sf::VideoMode(1280, 720, 32), "FreeGravity by Ivan Maklyakov (C) 2014",
                            sf::Style::Titlebar | sf::Style::Close, settings);

    double dt = 0.3;
    int fps = 0;
    int frames = 0;

    Gravity gravity(Gravity::INTEGRATION_RUNGE_KUTTA4, Gravity::OPTIMISATION_TREE);
    gravity.Add(2048);

    sf::Clock clock;
    float t0 = clock.getElapsedTime().asSeconds();

    sf::Text text("", font, 16);
    text.setPosition(4, 0);
    sf::Text about("Keys: Space - add particles, \"-\" / \"+\" - change dt, \"/\" / \"*\" - zoom, "
                   "Enter - reset view, Backspace - to center of mass",
                   font, 16);
    about.setPosition(4, 698);

    sf::View view = window.getView();
    sf::View view_default = window.getView();

    while (window.isOpen())
    {
        sf::Event event;
        while (window.pollEvent(event))
        {
            if (event.type == sf::Event::Closed)
                window.close();

            if (event.type == sf::Event::KeyPressed)
            {
                if (event.key.code == sf::Keyboard::Escape)
                    window.close();
                if (event.key.code == sf::Keyboard::Space && gravity.GetCount() < 4096)
                    gravity.Add(512);
                if (event.key.code == sf::Keyboard::Add && dt < 4.99)
                    dt += dt * 0.1;
                if (event.key.code == sf::Keyboard::Subtract && dt > 0.01)
                    dt -= dt * 0.1;

                if (event.key.code == sf::Keyboard::Multiply)
                    view.zoom(0.9f);
                if (event.key.code == sf::Keyboard::Divide)
                    view.zoom(1.1f);

                if (event.key.code == sf::Keyboard::BackSpace)
                {
                    float x, y;
                    gravity.GetCenter(x, y);
                    view.setCenter(x, y);
                }

                if (event.key.code == sf::Keyboard::Return)
                    view = view_default;
            }
        }

        float t1 = clock.getElapsedTime().asSeconds();
        gravity.DoGravity(dt);
        float t2 = clock.getElapsedTime().asSeconds();
        gravity.DoCollisions();
        float t3 = clock.getElapsedTime().asSeconds();
        gravity.DoClean();
        float t4 = clock.getElapsedTime().asSeconds();
        window.setView(view);
        window.clear(sf::Color(0, 0, 0, 255));
        gravity.Draw(window);
        frames++;
        float t5 = clock.getElapsedTime().asSeconds();

        if (t5 - t0 > 1.0f)
        {
            t0 += 1.0;
            fps = frames;
            frames = 0;
        }

        std::stringstream txt;
        txt << "FPS = " << fps << "\n";
        txt << "Count = " << gravity.GetCount() << "\n";
        txt << "dt = " << dt << "\n";
        txt << "DoGravity = " << (t2 - t1) * 1000.0f << "\n";
        txt << "DoCollisions = " << (t3 - t2) * 1000.0f << "\n";
        txt << "DoClean = " << (t4 - t3) * 1000.0f << "\n";
        txt << "Draw = " << (t5 - t4) * 1000.0f;
        text.setString(txt.str());

        window.setView(view_default);
        window.draw(text);
        window.draw(about);

        window.display();
    }

    return EXIT_SUCCESS;
}
