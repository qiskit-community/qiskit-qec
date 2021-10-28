#include <iostream>
#include <string>

class Entity
{
private:
    std::string m_Name;
    int m_Age;

public:
    Entity(const std::string &name) : m_Name(name), m_Age(-1) {}

    Entity(int age) : m_Name("Unknown"), m_Age(age) {}

    std::string const &GetName() const
    {
        return m_Name;
    }
};

void PrintEntity(const Entity &entity)
{
    std::cout << entity.GetName() << std::endl;
}

int main()
{
    Entity e0("cheeseburger");
    Entity e1(10);

    Entity e3 = std::string("Chernie");

    const std::string name = "Cigarettes";

    PrintEntity(name);
}