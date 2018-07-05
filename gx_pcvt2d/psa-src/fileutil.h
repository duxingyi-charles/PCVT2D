/**
 * This file is part of the point set analysis tool psa
 * 
 * Copyright 2011, Thomas Schl\"{o}mer, thomas.schloemer@uni-konstanz.de
 *
 * Thanks to Daniel Heck.
 * 
 * psa is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef FILEUTIL_H
#define FILEUTIL_H

#include <sys/types.h>
#include <sys/stat.h>
#include <string>

#if defined (_MSC_VER)

#ifndef S_ISREG
#define S_ISREG(m) ((m & _S_IFREG) == _S_IFREG)
#endif

#ifndef S_ISDIR
#define S_ISDIR(m) (((m) & _S_IFDIR) == _S_IFDIR)
#endif

#endif


/**
 *  Determines if the specified name exists and represents a file.
 */
inline bool FileExists(const std::string &fname) {
    struct stat s;
    return (stat(fname.c_str(), &s) == 0 && S_ISREG(s.st_mode));
}


/**
 *  Determines if the specified name exists and represents a folder.
 */
inline bool FolderExists(const std::string &fname) {
    struct stat s;
    return (stat(fname.c_str(), &s) == 0 && S_ISDIR(s.st_mode));
}


/**
 *  Remove trailing slash from a string if present
 */
inline void RemoveTrailingSlash(std::string &fname) {
#if defined (_MSC_VER)
    if ( fname[fname.length()-1] == '\\' || fname[fname.length()-1] == '/' ) {
        fname.erase(fname.length()-1, 1);
    }
#else
    if ( fname[fname.length()-1] == '/') {
        fname.erase(fname.length()-1, 1);
    }
#endif
}


/**
 *  Extract file name from a string
 */
inline std::string GetFileName(std::string &fname)
{
    size_t found = fname.find_last_of("/\\");

    if (found == std::string::npos)
        return fname;
    else
        return fname.substr(found + 1);
}


/**
 *  Directory iterators. Thanks to Daniel Heck.
 */
struct DirEntry
{
    std::string name;
    bool is_dir;
};

class DirIter
{
public:
    virtual ~DirIter () {}
    
    virtual bool open (const std::string &path) = 0;
    
    virtual bool get_next (DirEntry &entry) = 0;
};


/** -------------------- DirIter (Win32) -------------------- */

#if defined (_MSC_VER)

#include <windows.h>

class DirIterOS : DirIter
{
public:
    
    DirIterOS (const std::string &path) 
    : m_handle (INVALID_HANDLE_VALUE)
    {
        open (path);
    }
    
    ~DirIterOS () {
        close();
    }
    
    bool open (const std::string &path)
    {
        std::string glob (path);
        glob += "\\*.*";
        m_handle = FindFirstFile (glob.c_str(), &m_dir);
        return m_handle != INVALID_HANDLE_VALUE;
    }
    
    bool get_next (DirEntry &entry) {
        if (m_handle != INVALID_HANDLE_VALUE)
        {
            entry.name = m_dir.cFileName;
            entry.is_dir = m_dir.dwFileAttributes & FILE_ATTRIBUTE_DIRECTORY;
            if (!FindNextFile (m_handle, &m_dir)) 
                close();
            return true;
        }
        return false;
    }
    
private:
    
    void close ()
    {
        if (m_handle != INVALID_HANDLE_VALUE)
        {
            FindClose (m_handle);
            m_handle = INVALID_HANDLE_VALUE;
        }
    }
    
    // Variables.
    WIN32_FIND_DATA m_dir;
    HANDLE            m_handle;
};


/** -------------------- DirIter (POSIX) -------------------- */

#else

#include <dirent.h>

class DirIterOS : DirIter
{
public:
    
    DirIterOS (const std::string &path) : m_dir (NULL), m_entry (NULL) {
        open (path);
        dir_path = path;
    }
    
    ~DirIterOS () {
        if (m_dir != NULL)
            closedir (m_dir);
    }
    
    bool open (const std::string &path)
    {
        m_dir = opendir (path.c_str());
        return m_dir != 0;
    }
    
    bool get_next (DirEntry &entry)
    {
        if (m_dir == 0) return false;
        m_entry = readdir(m_dir);
        if (m_entry != NULL)
        {
            entry.name = m_entry->d_name;
            std::string fname = dir_path + "/" + entry.name;
            struct stat s;
            entry.is_dir = stat(fname.c_str(), &s)==0 && S_ISDIR(s.st_mode);
            return true;
        }
        return false;
    }
    
private:
    
    std::string    dir_path;
    DIR              *m_dir;
    struct dirent *m_entry;
};

#endif


#endif // FILEUTIL_H

