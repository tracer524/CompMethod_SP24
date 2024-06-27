#include <iostream>

class Point{
    public:
        double x, y;
};

class Vector{
    public:
        double x, y;
        void init_2(Point A, Point B);
        void init_1(Point A);
};

void Vector::init_2(Point A, Point B){
    x = B.x - A.x;
    y = B.y - A.y;
}

void Vector::init_1(Point A){
    x = A.x
}

double distance(Point A, Point B){
    double delta_x = A.x - B.x;
    double delta_y = A.y - B.y;
    return sqrt(delta_x * delta_x + delta_y * delta_y);
}

double angle(Vector A, Vector B){

}

int main(){
    return 0;
}