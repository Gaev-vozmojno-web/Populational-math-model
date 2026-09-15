import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import numpy as np


def a(i, tau):
    """Определение коэффициентов"""
    return 0.03 * float(i) * tau


def main(T, L, N, M, FI, C, U_max):
    t1 = T / M
    tau = L / N
    gamma = t1 / tau

    # Выводим гамму на экран для проверки устойчивости
    print(f"\n[Инфо] Текущее значение gamma = {gamma:.4f}")
    if gamma > 1.0:
        print("[ВНИМАНИЕ] gamma > 1! График развалится на иголки. Увеличьте M!")

    # Используем матрицу NumPy (размер N+1 x M+1)
    y = np.zeros((N + 1, M + 1))

    # Задаем начальные условия для первого столбца (j = 0) по всему пространству
    for i in range(N + 1):
        y[i, 0] = FI

    # Численно решаем
    for j in range(M):
        s = 0.0
        for i in range(N, 0, -1):  
            y[i, j + 1] = (
                (1 - gamma) * y[i, j]
                - a(i, tau) * t1 * U_max
                + gamma * y[i - 1, j]
            )

            if y[i, j + 1] < 0:
                y[i, j + 1] = 0.0

            # решаем интеграл
            if i == 0 or i == N:
                s = s + C * y[i, j] * tau / 2.0
            else:
                s = s + C * y[i, j] * tau

        # Расчет плотности
        y[0, j + 1] = 2.0 / (2.0 - tau * C) * s
        if y[0, j + 1] < 0:
            y[0, j + 1] = 0.0

    # 1. Создание матрицы и СОХРАНЕНИЕ результатов в текстовый файл
    Y_transposed = y.T  
    txt_filename = "result.txt"

    with open(txt_filename, "w") as f1:
        for j in range(M + 1):
            line = " ".join(f"{Y_transposed[j, i]:6.2f}" for i in range(N + 1))
            f1.write(line + "\n")

    # 2. Построение и сохранение 3D-графика под матрицу
    img_filename = "3d_graph.png"

    x_space = np.linspace(0, L, N + 1)  
    t_time = np.linspace(0, T, M + 1)  

    X_mesh, T_mesh = np.meshgrid(x_space, t_time)

    fig = plt.figure(figsize=(10, 7))
    ax = fig.add_subplot(111, projection="3d")

    # Включаем легкую темную сетку поверх графика, чтобы лучше видеть изгибы ("волны")
    surf = ax.plot_surface(
        X_mesh, T_mesh, Y_transposed, cmap="viridis", edgecolor="black", linewidth=0.1, alpha=0.9
    )

    ax.set_title("3D Визуализация популяционной модели $U(\\tau, t)$", fontsize=14, pad=15)
    ax.set_xlabel("Возраст ($\\tau$)", fontsize=12, labelpad=10)
    ax.set_ylabel("Время ($t$)", fontsize=12, labelpad=10)
    ax.set_zlabel("Плотность $U(\\tau, t)$", fontsize=12, labelpad=10)

    fig.colorbar(surf, ax=ax, shrink=0.5, aspect=10, label="$U(\\tau, t)$")
    ax.view_init(elev=25, azim=67)

    plt.tight_layout()
    plt.savefig(img_filename, dpi=300)
    plt.close()

    print(f"[Успешно] Расчет завершен!")
    print(f" -> Картинка 3D-графика: {img_filename}")


if __name__ == "__main__":
    user_T = float(input("Введите T: "))
    user_L = float(input("Введите L: "))
    user_N = int(input("Введите N: "))
    user_M = int(input("Введите M: "))
    user_FI = float(input("Введите фи: "))
    user_C = float(input("Введите с: "))
    user_U_max = float(input("Введите Umax: "))
    
    main(user_T, user_L, user_N, user_M, user_FI, user_C, user_U_max)




