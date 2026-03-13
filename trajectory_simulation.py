"""
Численное моделирование смещения траектории ракеты при запуске с морской платформы
===================================================================================

Задача: Оцените смещение траектории на высоте 10 км из-за качки морской платформы
        (период 8 с, амплитуда 2°, задержка управления 0.3 с).

Метод: Численное интегрирование уравнений движения методом Рунге-Кутты 4-го порядка
       (Метод B из списка литературы)

Источники:
  [1] Лебедев А.А., Чернобровкин Л.С. — Динамика полёта БПЛА — М.: Машиностроение, 1973
      (уравнения движения, гл. 2–3)
  [2] Колесников К.С. — Динамика ракет — М.: Машиностроение, 2003
      (начальные условия от качки подвижного основания)
  [3] Blakelock J.H. — Automatic Control of Aircraft and Missiles — Wiley, 1991
      (задержка в контуре управления, передаточные функции автопилота)
  [7] Frosch J.A., Vallely D.P. — Saturn AS-501/S-IC Flight Control System Design
      Journal of Spacecraft and Rockets, Vol. 4, No. 8, 1967
      (методология численного моделирования, TVC, шаг интегрирования)

Автор: Численная модель «запуск с морской платформы»
Версия: 1.0
"""

import os
import matplotlib
matplotlib.use('Agg')   # headless-режим: установить до импорта pyplot
import numpy as np
import matplotlib.pyplot as plt
from collections import deque


# =============================================================================
# ПАРАМЕТРЫ ЗАДАЧИ (источник: условие задачи + [2] Колесников)
# =============================================================================
T_PLATFORM = 8.0              # с,   период качки платформы
A_PLATFORM = np.radians(2.0)  # рад, амплитуда качки платформы (2°)
TAU_DELAY  = 0.3              # с,   задержка в контуре управления

# =============================================================================
# ПАРАМЕТРЫ РАКЕТЫ (типичная орбитальная РН среднего класса — Зенит/Sea Launch)
# источник: [4] Sutton, открытые данные Sea Launch
# =============================================================================
M0        = 459_600.0  # кг, начальная масса ракеты
M_PROPEL  = 362_000.0  # кг, масса топлива первой ступени
P_THRUST  = 7_257_000  # Н,  тяга двигателя первой ступени (РД-171)
ISP       = 309.0      # с,  удельный импульс (РД-171, атмосферный)
G0        = 9.80665    # м/с², стандартное ускорение свободного падения
M_DOT     = P_THRUST / (ISP * G0)  # ~2396 кг/с, расход топлива

L_ROCKET  = 59.6       # м,  длина первой ступени
D_ROCKET  = 3.9        # м,  диаметр
S_REF     = np.pi * (D_ROCKET / 2) ** 2  # м², площадь миделя (~11.9 м²)

# Аэродинамика (приближённые значения)
CN_ALPHA = 10.0    # 1/рад, производная нормальной силы по углу атаки
CD0      = 0.25    # безразмерный, коэффициент лобового сопротивления (дозвук)

# Геометрия (положения от носа ракеты, положительное — к хвосту)
XCG0   = 0.42 * L_ROCKET  # м, начальное положение ЦМ
XCP    = 0.55 * L_ROCKET  # м, положение ЦД (аэродинамически неустойчива: XCP > XCG)
L_TVC  = 0.50 * L_ROCKET  # м, плечо TVC (от ЦМ до шарнира сопла)

# Момент инерции (приближение: цилиндр с неравномерной заправкой)
def moment_of_inertia(mass):
    return mass * L_ROCKET ** 2 / 14.0  # кг·м², коэффициент ~1/14 для ракеты с баком

# Параметры TVC (источник: [4] Sutton, [7] Frosch)
DELTA_MAX = np.radians(6.0)   # рад, максимальное отклонение сопла (6°)

# =============================================================================
# ПАРАМЕТРЫ АВТОПИЛОТА (источник: [3] Blakelock, [7] Frosch)
# =============================================================================
KP = 2.5   # коэффициент пропорциональной части (рад/рад)
KD = 0.8   # коэффициент дифференциальной части (с/рад)


# =============================================================================
# МОДЕЛЬ АТМОСФЕРЫ (стандартная ISA-1976, экспоненциальное приближение)
# =============================================================================
def atmosphere(h_m):
    """Плотность воздуха [кг/м³] на высоте h_m [м] — экспоненциальное приближение ISA."""
    rho0   = 1.225    # кг/м³, плотность у земли
    H_SCALE = 8_435.0  # м,    масштабная высота
    return rho0 * np.exp(-h_m / H_SCALE)


# =============================================================================
# ДВИЖЕНИЕ ПЛАТФОРМЫ (источник: условие задачи + [2] Колесников)
# =============================================================================
OMEGA_PLAT = 2.0 * np.pi / T_PLATFORM  # рад/с, угловая частота качки


def platform_angle(t, phase=0.0):
    """Угол платформы в момент времени t [с], начальная фаза phase [рад]."""
    return A_PLATFORM * np.cos(OMEGA_PLAT * t + phase)


def platform_rate(t, phase=0.0):
    """Угловая скорость платформы [рад/с]."""
    return -A_PLATFORM * OMEGA_PLAT * np.sin(OMEGA_PLAT * t + phase)


# =============================================================================
# УРАВНЕНИЯ ДВИЖЕНИЯ (источник: [1] Лебедев-Чернобровкин, гл. 2–3)
# =============================================================================
def equations_of_motion(t, state, delta_nozzle, mass):
    """
    Уравнения движения в плоскости «бок-высота» (упрощённые 3DOF).

    Вектор состояния state = [y, vy, h, vh, phi, omega]:
        y      — боковое смещение [м]
        vy     — боковая скорость [м/с]
        h      — высота [м]
        vh     — вертикальная скорость [м/с]
        phi    — угол тангажа от вертикали (+ в сторону +y) [рад]
        omega  — угловая скорость тангажа [рад/с]

    delta_nozzle — угол отклонения сопла (+ отклоняет в сторону +y) [рад]
    mass         — текущая масса [кг]

    Возвращает: d(state)/dt
    """
    y, vy, h, vh, phi, omega = state

    # Ограничение отклонения сопла
    delta = np.clip(delta_nozzle, -DELTA_MAX, DELTA_MAX)

    # Скорость и динамическое давление
    V = np.sqrt(vy ** 2 + vh ** 2)
    if V < 0.1:
        V = 0.1
    rho = atmosphere(h)
    q = 0.5 * rho * V ** 2

    # Угол атаки (угол между вектором скорости и продольной осью ракеты)
    gamma = np.arctan2(vy, vh)   # угол траектории от вертикали
    alpha = phi - gamma           # угол атаки

    # -----------------------------------------------------------------------
    # Силы
    # -----------------------------------------------------------------------
    # Тяга: направлена вдоль оси ракеты, сопло отклонено на delta
    # В инерциальных осях (y, h):
    #   Вектор оси ракеты: (sin(phi), cos(phi))
    #   Тяга со смещением сопла: поворот на -delta от оси
    Thrust_h = P_THRUST * np.cos(phi - delta)
    Thrust_y = P_THRUST * np.sin(phi - delta)

    # Нормальная аэродинамическая сила (перпендикулярна оси ракеты)
    F_N = CN_ALPHA * alpha * q * S_REF
    # В инерциальных осях:
    FN_h = -F_N * np.sin(phi)
    FN_y =  F_N * np.cos(phi)

    # Аэродинамическое сопротивление (вдоль вектора скорости, против движения)
    drag = CD0 * q * S_REF
    Drag_h = -drag * (vh / V)
    Drag_y = -drag * (vy / V)

    # -----------------------------------------------------------------------
    # Ускорения ЦМ
    # -----------------------------------------------------------------------
    ay = (Thrust_y + FN_y + Drag_y) / mass
    ah = (Thrust_h + FN_h + Drag_h) / mass - G0

    # -----------------------------------------------------------------------
    # Угловое ускорение (вокруг ЦМ)
    # -----------------------------------------------------------------------
    Izz = moment_of_inertia(mass)

    # Момент от TVC (на плече L_TVC; положительный delta создаёт отрицательный момент)
    M_tvc  = -P_THRUST * L_TVC * np.sin(delta)

    # Момент аэродинамических сил (по Лебедеву: нестабилизирующий для ракеты)
    arm    = XCG0 - XCP   # < 0 → ракета аэродинамически неустойчива
    M_aero = F_N * arm

    alpha_dot = (M_tvc + M_aero) / Izz

    return [vy, ay, vh, ah, omega, alpha_dot]


# =============================================================================
# БУФЕР ЗАДЕРЖКИ (источник: [3] Blakelock; аппаратная реализация τ = 0.3 с)
# =============================================================================
class DelayBuffer:
    """
    Дискретный буфер задержки для имитации запаздывания в контуре управления.
    Хранит историю значений за последние tau секунд с шагом dt.
    """
    def __init__(self, tau, dt, initial_value=0.0):
        self.tau = tau
        self.dt  = dt
        n_steps  = max(1, int(round(tau / dt)))
        self.buf = deque([initial_value] * (n_steps + 1), maxlen=(n_steps + 1))
        self.n_steps = n_steps

    def push(self, value):
        """Добавить новое значение."""
        self.buf.append(value)

    def get_delayed(self):
        """Вернуть значение tau секунд назад."""
        return self.buf[0]


# =============================================================================
# RK4 ИНТЕГРАТОР (источник: [7] Frosch, шаг 0.01 с)
# =============================================================================
def rk4_step(f, t, state, dt, *args):
    """Один шаг метода Рунге-Кутты 4-го порядка."""
    k1 = np.array(f(t,          state,              *args))
    k2 = np.array(f(t + dt/2,   state + dt/2 * k1,  *args))
    k3 = np.array(f(t + dt/2,   state + dt/2 * k2,  *args))
    k4 = np.array(f(t + dt,     state + dt  * k3,   *args))
    return state + (dt / 6.0) * (k1 + 2*k2 + 2*k3 + k4)


# =============================================================================
# ОСНОВНАЯ ПРОЦЕДУРА МОДЕЛИРОВАНИЯ
# =============================================================================
def run_simulation(scenario='max_angle', dt=0.01, h_target=10_000.0):
    """
    Запустить моделирование траектории.

    Parameters
    ----------
    scenario : str
        'max_angle' — старт при максимальном угле платформы (φ₀ = 2°, ω₀ = 0)
        'max_rate'  — старт при максимальной угловой скорости платформы (φ₀ = 0, ω₀ = max)
    dt : float
        Шаг интегрирования [с]. По умолчанию 0.01 с (как у Frosch [7]).
    h_target : float
        Высота, на которой снимается результат [м]. По умолчанию 10 000 м.

    Returns
    -------
    dict с результатами: временные ряды и ключевые результаты
    """
    # --- Начальные условия из динамики платформы ([2] Колесников) ---
    if scenario == 'max_angle':
        # Платформа в момент максимального отклонения: φ₀=2°, ω₀=0
        phase   = 0.0
        phi_0   = A_PLATFORM                # рад
        omega_0 = 0.0                       # рад/с
    else:  # 'max_rate'
        # Платформа при максимальной угловой скорости поворота: φ₀=0, ω₀=max
        phase   = np.pi / 2.0
        phi_0   = 0.0                        # рад
        omega_0 = -A_PLATFORM * OMEGA_PLAT   # рад/с (максимальная скорость)

    # --- Начальное состояние [y, vy, h, vh, phi, omega] ---
    state = np.array([0.0, 0.0, 0.0, 0.0, phi_0, omega_0])
    mass  = M0

    # --- Буфер задержки для φ и ω (источник: [3] Blakelock) ---
    delay_phi   = DelayBuffer(TAU_DELAY, dt, initial_value=phi_0)
    delay_omega = DelayBuffer(TAU_DELAY, dt, initial_value=omega_0)

    # --- Хранение истории ---
    t = 0.0
    t_arr     = [t]
    y_arr     = [state[0]]
    h_arr     = [state[2]]
    phi_arr   = [state[4]]
    delta_arr = [0.0]
    mass_arr  = [mass]

    y_at_target = None
    t_at_target = None

    # --- Цикл интегрирования ---
    while True:
        # Обновить буфер задержки
        delay_phi.push(state[4])
        delay_omega.push(state[5])

        # Получить задержанные измерения (именно то, что «видит» автопилот)
        phi_delayed   = delay_phi.get_delayed()
        omega_delayed = delay_omega.get_delayed()

        # Закон управления (PD-регулятор; источник: [3] Blakelock, [7] Frosch)
        # Цель: удерживать ракету вертикально (phi_ref = 0)
        # Знак: при phi > 0 нужен delta > 0, чтобы момент M_tvc = -P*L*sin(delta) < 0
        # уменьшал phi, и тяга Thrust_y = P*sin(phi-delta) тянула в -y
        delta_cmd = (KP * phi_delayed + KD * omega_delayed)
        delta_cmd = np.clip(delta_cmd, -DELTA_MAX, DELTA_MAX)

        # Шаг RK4
        state = rk4_step(equations_of_motion, t, state, dt, delta_cmd, mass)

        # Обновить массу (уменьшается с постоянной скоростью)
        mass -= M_DOT * dt
        m_dry = M0 - M_PROPEL
        if mass < m_dry:
            mass = m_dry   # Топливо закончилось

        t += dt

        # Защита от ухода за пределы
        if state[2] < -10:   # Ракета упала
            break

        # Проверка достижения целевой высоты
        if state[2] >= h_target and y_at_target is None:
            y_at_target = float(state[0])
            t_at_target = t
            # Продолжать симуляцию нет смысла
            break

        # Запись истории
        t_arr.append(t)
        y_arr.append(state[0])
        h_arr.append(state[2])
        phi_arr.append(state[4])
        delta_arr.append(delta_cmd)
        mass_arr.append(mass)

        # Максимальное время симуляции (предохранитель)
        if t > 200.0:
            break

    return {
        't':            np.array(t_arr),
        'y':            np.array(y_arr),
        'h':            np.array(h_arr),
        'phi':          np.array(phi_arr),
        'delta':        np.array(delta_arr),
        'mass':         np.array(mass_arr),
        'y_at_target':  y_at_target,
        't_at_target':  t_at_target,
        'scenario':     scenario,
    }


# =============================================================================
# ЧАСТОТНАЯ ЭКСПРЕСС-ОЦЕНКА (Метод C, источник: [3] Blakelock, [8] Smith)
# =============================================================================
def frequency_estimate():
    """
    Быстрая аналитическая оценка фазового сдвига и порядка ошибки.
    Источник: [3] Blakelock (передаточные функции автопилота),
              [8] Smith (предиктор и анализ задержки).
    """
    omega = 2.0 * np.pi / T_PLATFORM   # рад/с
    tau   = TAU_DELAY                  # с
    phi_shift = omega * tau            # рад
    phi_shift_deg = np.degrees(phi_shift)

    # Угловое отклонение в момент реакции (при максимальной скорости качки)
    delta_phi_max = A_PLATFORM * omega * tau  # рад

    # Оценка начальной поперечной скорости (через время задержки * ускорение)
    # a_lateral = P * phi_0 / m (при phi_0 = delta_phi_max)
    a_lat = P_THRUST * delta_phi_max / M0
    v_lat = a_lat * tau   # м/с, приближённая начальная поперечная скорость

    return {
        'phase_shift_rad': phi_shift,
        'phase_shift_deg': phi_shift_deg,
        'delta_phi_max_deg': np.degrees(delta_phi_max),
        'initial_lateral_accel': a_lat,
        'initial_lateral_speed': v_lat,
    }


# =============================================================================
# ВЫВОД РЕЗУЛЬТАТОВ И ПОСТРОЕНИЕ ГРАФИКОВ
# =============================================================================
def print_results(res_max_angle, res_max_rate):
    """Вывести таблицу результатов в консоль."""
    print()
    print("=" * 65)
    print("  РЕЗУЛЬТАТЫ МОДЕЛИРОВАНИЯ: Запуск с морской платформы")
    print("=" * 65)
    print(f"  Платформа: T = {T_PLATFORM} с, A = {np.degrees(A_PLATFORM):.1f}°")
    print(f"  Задержка управления: τ = {TAU_DELAY} с")
    print(f"  Целевая высота: 10 000 м")
    print()

    fe = frequency_estimate()
    print("  [Экспресс-оценка, Метод C — источник: Blakelock [3], Smith [8]]")
    print(f"  Фазовый сдвиг ωτ = {fe['phase_shift_deg']:.2f}°")
    print(f"  Непоправленное угловое отклонение Δφ ≈ {fe['delta_phi_max_deg']:.3f}°")
    print()

    print("  [Точный расчёт, Метод B — RK4, шаг dt = 0.01 с]")
    print(f"  {'Сценарий':<30} {'t до 10 км, с':>13} {'Боковое смещение, м':>20}")
    print("  " + "-" * 63)

    for res in [res_max_angle, res_max_rate]:
        scen_name = {
            'max_angle': 'Макс. угол (φ₀=2°, ω₀=0)',
            'max_rate':  'Макс. скорость (φ₀=0°, ω₀=max)',
        }[res['scenario']]

        if res['y_at_target'] is not None:
            y_str = f"{res['y_at_target']:+.1f}"
            t_str = f"{res['t_at_target']:.1f}"
        else:
            y_str = "не достигнута"
            t_str = "—"

        print(f"  {scen_name:<30} {t_str:>13} {y_str:>20}")

    print()

    # Наихудший случай
    valid_results = [r for r in [res_max_angle, res_max_rate] if r['y_at_target'] is not None]
    if valid_results:
        worst = max(valid_results, key=lambda r: abs(r['y_at_target']))
        y_worst = worst['y_at_target']

        print("  НАИХУДШИЙ СЛУЧАЙ (детерминированная верхняя граница):")
        print(f"  Максимальное боковое смещение на h = 10 км: {abs(y_worst):.0f} м")
        print()
        print("  ВЫВОД:")
        if abs(y_worst) < 1000:
            print(f"  Смещение {abs(y_worst):.0f} м << 1000 м (допуск ±1 км).")
            print("  ✅ Смещением МОЖНО ПРЕНЕБРЕЧЬ для вывода спутника на орбиту.")
        else:
            print(f"  Смещение {abs(y_worst):.0f} м > 1000 м (допуск ±1 км).")
            print("  ⚠️  Смещение НЕЛЬЗЯ ИГНОРИРОВАТЬ — требуется коррекция траектории.")

    print("=" * 65)
    print()


def plot_results(res_max_angle, res_max_rate, save_path=None):
    """Построить графики траектории."""
    fig, axes = plt.subplots(2, 2, figsize=(12, 9))
    fig.suptitle(
        f"Смещение траектории с морской платформы\n"
        f"(T={T_PLATFORM}с, A=2°, τ={TAU_DELAY}с, Метод: RK4)",
        fontsize=13
    )

    colors = {'max_angle': 'steelblue', 'max_rate': 'darkorange'}
    labels = {
        'max_angle': r'$\varphi_0 = 2°,\; \dot\varphi_0 = 0$',
        'max_rate':  r'$\varphi_0 = 0°,\; \dot\varphi_0 = \max$',
    }

    for res in [res_max_angle, res_max_rate]:
        c = colors[res['scenario']]
        lbl = labels[res['scenario']]
        t, y, h, phi, delta = res['t'], res['y'], res['h'], res['phi'], res['delta']

        # 1) Траектория y(h)
        axes[0, 0].plot(y, h / 1e3, color=c, label=lbl)
        axes[0, 0].set_xlabel('Боковое смещение y [м]')
        axes[0, 0].set_ylabel('Высота h [км]')
        axes[0, 0].set_title('Траектория y(h)')
        axes[0, 0].axhline(10, ls='--', color='gray', lw=0.8, label='h = 10 км')
        axes[0, 0].legend(fontsize=9)
        axes[0, 0].grid(True, alpha=0.4)

        # 2) Боковое смещение y(t)
        axes[0, 1].plot(t, y, color=c, label=lbl)
        axes[0, 1].set_xlabel('Время [с]')
        axes[0, 1].set_ylabel('Боковое смещение y [м]')
        axes[0, 1].set_title('Смещение y(t)')
        axes[0, 1].legend(fontsize=9)
        axes[0, 1].grid(True, alpha=0.4)

        # 3) Угол тангажа φ(t)
        axes[1, 0].plot(t, np.degrees(phi), color=c, label=lbl)
        axes[1, 0].set_xlabel('Время [с]')
        axes[1, 0].set_ylabel('Угол тангажа φ [°]')
        axes[1, 0].set_title('Угол тангажа φ(t)')
        axes[1, 0].axhline(0, ls='--', color='gray', lw=0.8)
        axes[1, 0].legend(fontsize=9)
        axes[1, 0].grid(True, alpha=0.4)

        # 4) Отклонение сопла δ(t)
        axes[1, 1].plot(t, np.degrees(delta), color=c, label=lbl)
        axes[1, 1].set_xlabel('Время [с]')
        axes[1, 1].set_ylabel('Отклонение сопла δ [°]')
        axes[1, 1].set_title('Команда TVC δ(t)')
        axes[1, 1].axhline(6,   ls=':', color='red', lw=0.8, label='Ограничение ±6°')
        axes[1, 1].axhline(-6,  ls=':', color='red', lw=0.8)
        axes[1, 1].legend(fontsize=9)
        axes[1, 1].grid(True, alpha=0.4)

    plt.tight_layout()
    if save_path:
        plt.savefig(save_path, dpi=150, bbox_inches='tight')
        print(f"  График сохранён: {save_path}")
    else:
        plt.show()


# =============================================================================
# ТОЧКА ВХОДА
# =============================================================================
if __name__ == '__main__':
    print("Запуск моделирования... (Метод RK4, dt=0.01 с)")
    print("Источники: [1] Лебедев, [2] Колесников, [3] Blakelock, [7] Frosch")

    # Запустить оба сценария (наихудший угол и наихудшая скорость)
    res_angle = run_simulation(scenario='max_angle', dt=0.01)
    res_rate  = run_simulation(scenario='max_rate',  dt=0.01)

    # Вывести результаты в консоль
    print_results(res_angle, res_rate)

    # Сохранить график
    out_dir = os.path.dirname(os.path.abspath(__file__))
    plot_path = os.path.join(out_dir, 'trajectory_plot.png')
    try:
        plot_results(res_angle, res_rate, save_path=plot_path)
    except Exception as e:
        print(f"  (График не построен: {e})")
