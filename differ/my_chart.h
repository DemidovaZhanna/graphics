#ifndef MY_CHART_H
#define MY_CHART_H

#include <QtCharts>

class MyChart : public QChart {
	Q_OBJECT
public:
	MyChart(int a, int b);
	void addSeries(QAbstractSeries *a);

private:
	QValueAxis *axisX, *axisY;
};



#endif // MY_CHART_H
