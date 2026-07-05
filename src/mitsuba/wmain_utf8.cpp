/*
    This file is part of Mitsuba, a physically based rendering system.

    Copyright (c) 2007-2014 by Wenzel Jakob and others.

    Mitsuba is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License Version 3
    as published by the Free Software Foundation.

    Mitsuba is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program. If not, see <http://www.gnu.org/licenses/>.
*/

#if defined(_WIN32)

#include <Windows.h>
#include <vector>

extern int mts_main(int argc, char **argv);

namespace {

class UTF8CommandLine {
public:
    UTF8CommandLine(int argc, wchar_t *wargv[]) : m_argc(-1), m_argv(NULL), m_data(NULL) {
        if (argc <= 0)
            return;

        m_argc = argc;
        m_argv = new char *[argc];
        std::vector<int> lengths(argc);
        int total = 0;

        for (int i = 0; i < argc; ++i) {
            const int length = WideCharToMultiByte(CP_UTF8, 0, wargv[i], -1,
                                                   NULL, 0, NULL, NULL);
            if (length == 0) {
                m_argc = i;
                return;
            }

            total += length;
            lengths[i] = length;
        }

        m_data = new char[total];
        int offset = 0;

        for (int i = 0; i < m_argc; ++i) {
            int length = lengths[i];
            m_argv[i] = m_data + offset;
            length = WideCharToMultiByte(CP_UTF8, 0, wargv[i], -1,
                                         m_argv[i], length, NULL, NULL);
            if (length == 0) {
                m_argc = i;
                return;
            }

            offset += length;
        }
    }

    ~UTF8CommandLine() {
        delete[] m_argv;
        delete[] m_data;
    }

    int argc() const { return m_argc; }
    char **argv() const { return m_argv; }

private:
    int m_argc;
    char **m_argv;
    char *m_data;
};

} // namespace

int wmain(int argc, wchar_t *wargv[], wchar_t *) {
    UTF8CommandLine commandLine(argc, wargv);
    return mts_main(commandLine.argc(), commandLine.argv());
}

#endif
